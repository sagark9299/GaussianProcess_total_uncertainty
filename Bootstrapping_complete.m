clear
clc
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use bootstrapping to determine FIR surface uncertainty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Source GP module & FIR data
addpath('../GP_module/')
addpath('../New_FIR_data')
addpath('../Identification')
addpath('../Identification/data_time_series')
addpath('../Identification/TFDtools')
%% Identify flame model and load FIR data points
no_FIR_coeff = 25; %22
FIR_sampling_time = 4.16e-4; % 3.5e-4
no_temp_points = 13;
FIR_data = identification_FIR(no_FIR_coeff,FIR_sampling_time,no_temp_points);
%% Load FIR data
temp_space = [700,720,740,760,786,800,823,842,860,885,905,926,946]; %700,960
temperature_vec = temp_space(1:no_temp_points);
FIR_data = FIR_data(1:no_temp_points);
delta = FIR_data{1}.Ts;
time = 3*delta:delta:(no_FIR_coeff+3)*delta;

%% Splitting data into training and testing
training_index = [1,2,3,4,5,6,7,8,9,10,11,12,13];
% training_index = [1,2,3,4,5,6,7];
% testing_temperature = [823];
testing_index = [3]; 

temperature_vec(testing_index)
% Normalizing procedure
begin_time = time(1);
end_time = time(end);
time = (time-begin_time)/(end_time-begin_time);

temp_lower = min(temperature_vec(training_index));
temp_upper = max(temperature_vec(training_index));
temperature_vec = (temperature_vec-temp_lower)/(temp_upper-temp_lower);

%% Training data (nominal)
training_X = [];
training_y = [];

for i = 1:length(training_index)
    % Mean
    FIR = FIR_data{training_index(i)}.Numerator'; % nominal coefficient values
    % Construct training matrix
    training_X = [training_X;time',temperature_vec(training_index(i))*ones(length(time),1)];
    training_y = [training_y;FIR]; 
end
GP_FIR_nominal = GP(training_X,training_y,[],1,'off'); % [1.8698,1.4557] 1.698,1.1557

% Predictions
% pred_temp = temperature_vec(testing_index(1));
bootstrap = [100,100];    % 1:global realizations 2: local realizations
training_y_MC = zeros(bootstrap(1),size(training_y,1));
FIR_container = zeros(prod(bootstrap),length(time));
FIR_total = [];

pred_temp = lhsdesign(200,1);
tic

for k = 1:length(pred_temp)
[FIR_nominal,~] = pred_GP([time',pred_temp(k)*ones(length(time),1)],GP_FIR_nominal,'full');


%% Training data (realizations)
% bootstrap = [500,500];    % 1:global realizations 2: local realizations
% training_y_MC = zeros(bootstrap(1),size(training_y,1));
% FIR_container = zeros(prod(bootstrap),length(time));

% Generate realizations for training data
% FIR_total = [];
for i = 1:length(training_index)
    mu = FIR_data{training_index(i)}.Numerator';
    sigma = getcov(FIR_data{training_index(i)});
    training_y_MC(:,(i-1)*length(time)+1:i*length(time)) = mvnrnd(mu,sigma(1:end-1,1:end-1),bootstrap(1));
end

% Loop over each realization of the training data
for i = 1:bootstrap(1)
    % Obtain the updated GP model
    GP_FIR_MC = GP_NoTraining(training_X,training_y_MC(i,:)',GP_FIR_nominal.Theta);
    % Predict
    [FIR_MC,FIR_MC_cov] = pred_GP([time',pred_temp(k)*ones(length(time),1)],GP_FIR_MC,'full');
    FIR_MC_cov = (FIR_MC_cov+FIR_MC_cov')/2; % to ensure positive semi-definiteness
    % Realizations
    FIR_container((i-1)*bootstrap(2)+1:i*bootstrap(2),:) = mvnrnd(FIR_MC,FIR_MC_cov,bootstrap(2));
    
end
% temp = FIR_container;
FIR_total = [FIR_total;FIR_container];
end
toc
% Determine FIR statistics
FIR_nominal_uncertainty = cov(FIR_total);
FIR_pseudo_mean = mean(FIR_total);
%%
FIR_nominal_uncertainty = cov(FIR_total);
figure(10)
hold on
% Reference
% test_cov = getcov(FIR_data{testing_index(1)});
% p1 = plot(time*(end_time-begin_time)+begin_time,FIR_data{testing_index(1)}.Numerator,'o','MarkerFaceColor','r');
% plot(time*(end_time-begin_time)+begin_time,...
%     FIR_data{testing_index(1)}.Numerator'+1.96*sqrt(diag(test_cov(1:end-1,1:end-1))),'r--','LineWidth',2)
% plot(time*(end_time-begin_time)+begin_time,...
%     FIR_data{testing_index(1)}.Numerator'-1.96*sqrt(diag(test_cov(1:end-1,1:end-1))),'r--','LineWidth',2)

% GP-prediction
p2 = plot(time*(end_time-begin_time)+begin_time,FIR_pseudo_mean','-k','LineWidth',2);  
plot(time*(end_time-begin_time)+begin_time,FIR_pseudo_mean'+1.96*sqrt(diag(FIR_nominal_uncertainty)),'k--','LineWidth',2)  
plot(time*(end_time-begin_time)+begin_time,FIR_pseudo_mean'-1.96*sqrt(diag(FIR_nominal_uncertainty)),'k--','LineWidth',2)
% hold off
% legend([p1, p2],'Reference Case', 'GP Prediction');
% xlabel('Time (s)','FontSize',14);
% ylabel('Normalized FIR Coefficient','FontSize',14);
xlabel('$Time \ Lag, \  \tau, \ (s)$','FontSize',14,'Interpreter','Latex');
ylabel('$Normalized\ FIR\ Coefficient$','FontSize',14,'Interpreter','Latex');
ax = gca;
ax.FontSize = 12;
% axis([-0.3 0.5 0 0.009])
ax.XAxis.Exponent = 0;
xlim([0 0.01]);
%% Compare PDF of a single FIR coefficient
figure(2);
load ale.mat
load ale8.mat
load ale7.mat
load ale10.mat
% mu = 0.0642;
mu = FIR_pseudo_mean(15);
std_temp = FIR_total(:,15);
sigma = std(std_temp);
x = mu-(4*sigma):0.01:mu+(4*sigma);
y = normpdf(x,mu,sigma);
plot(x,y,'r',ale.x,ale.y,'--b',ale8.x8,ale8.y8,'--k',ale7.x7,ale7.y7,'--c',ale10.x10,ale10.y10,'--g','LineWidth',2)
xlabel('h (15)','FontSize',14);
ylabel('PDF','FontSize',14);
% legend('Total Uncertainty','Aleatoric Uncertainty');
ax = gca;
ax.FontSize = 12;