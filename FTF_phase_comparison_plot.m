clc
clear
close all hidden

%%
addpath('./data_time_series')
addpath('./TFDtools')
load experimental_FIR.mat % Experimental FIR obtained from FTF
% load('data_sagar_30pc_excitation_350ms_raw_rvel.mat') % base case i.e 823K
load('885K.mat');
data_base = data;
load('700K.mat');
data_960 = data;

load FTF_gain_10_new_arc.mat % Experimental data - Gain
load FTF_phase_10_new_arc.mat % Experimental data - Phase

%%
data_re_base = siResample(data_base,1/(2*1200)); % 1200 for FIR
data_norm_base = siNormalize(data_re_base);
data_cut_new_base = siCutSignal(data_norm_base,5e-3,350e-3);
rdata_base = resample(data_cut_new_base,1,1); % downsampling for cut-off frequency

%% Identification (FIR) or BJ
nb_FIR = 25; % number of non-zero impulse response coefficients % 
% 35 for base case
% 25 for best fit @ 860K
nk_FIR = 4; % number of time delays % 4 for base case
opt_FIR = impulseestOptions('RegularizationKernel','SS');
model_FIR = impulseest(rdata_base, nb_FIR,nk_FIR,opt_FIR); 

nb_BJ = 30; % 22 for current case 30; 35 for base case, 
% 30 for 860 K, where low frequency limit is 0.5
% 23 for 860 K, where low frequency limit is 1
nc_BJ = 10;% BJ parameters for error identification (numerator) 10 for base case
nd_BJ = 6; % BJ parameters for error identification (denominator) 10 for base case
nf_BJ = 0; % BJ paramter for deterministic identification
nk_BJ = 0;
opt_BJ = bjOptions('Display','off');
opt_BJ.Regularization.Lambda = 0.2;
R = eye(nb_BJ+nc_BJ+nd_BJ);
opt_BJ.Regularization.R = R;
model_BJ = bj(rdata_base,[nb_BJ nc_BJ nd_BJ nf_BJ nk_BJ],opt_BJ);

cov_FIR = getcov(model_FIR);
cov_BJ = getcov(model_BJ);

% model_container = cell(2,1);
% model_container{1} = model_FIR;
% model_container{2} = model_BJ;

delta = rdata_base.Ts;
% time = 0:delta:30*delta;
time = 3*delta:delta:(nb_FIR+3)*delta;
time_BJ = 0:delta:nb_BJ*delta;
% FIR
% stem(time(1:end-1),model_FIR.Numerator(2:end),'k','filled');
h1 = plot(time(1:nb_FIR),model_FIR.Numerator(2:end),'k','LineWidth',2);

hold on
plot(time(1:nb_FIR),model_FIR.Numerator(2:end)'+1.96*sqrt(diag(cov_FIR(2:end-1,2:end-1))),'--k','LineWidth',1.);
plot(time(1:nb_FIR),model_FIR.Numerator(2:end)'-1.96*sqrt(diag(cov_FIR(2:end-1,2:end-1))),'--k','LineWidth',1.);

% Box-Jenkins
% stem(time_BJ(1:end-1),model_BJ.B,'r','filled');
% h2 = plot(time_BJ(1:end-1),model_BJ.B,'r','LineWidth',2);

% plot(time_BJ(1:end-1),model_BJ.B'+1.96*sqrt(diag(cov_BJ(1:nb_BJ,1:nb_BJ))),'--r','LineWidth',1.);
% plot(time_BJ(1:end-1),model_BJ.B'-1.96*sqrt(diag(cov_BJ(1:nb_BJ,1:nb_BJ))),'--r','LineWidth',1.);

% Experiment
% stem(exp.time(1:30),exp.FIR(1:30),'c','filled')
h3 = plot(exp.time(1:nb_FIR),exp.FIR(1:nb_FIR),'c','LineWidth',2);
legend([h1 h3],'FIR','Experiment')
xlabel('Time (s)')
ylabel('Normalized FIR coefficient (-)')
ax = gca;
ax.FontSize = 12;
xlim([0 0.012])
xticks(0:0.002:0.012)
%% %--------Error bar from Experiments--------%
err_g = zeros(length(GainFDFxyscan(:,2)),1);
err_ph = zeros(length(PhaseFDFxyscan(:,2)),1);

err_g(1:29) = 0.1;
err_g(29:end) = linspace(0.1,0.3,length(err_g(29:end)));

err_ph(1:9) = 0.08*pi;
err_ph(10:13) = linspace(0.08*pi,0.3*pi,length(err_ph(10:13)));
err_ph(13:19) = linspace(0.3*pi,0.15*pi,length(err_ph(13:19)));
err_ph(20:end) = 0.15*pi;

% % %--------Gain--------%

[mag_FIR,phase_FIR,wout_FIR,sdmag_FIR,sdphase_FIR] = bode(model_FIR);
[mag_BJ,phase_BJ,wout_BJ,sdmag_BJ,sdphase_BJ] = bode(model_BJ);
C = colormap('lines');
err_c = C(1,:);
err_lw = 0.2;

figure(4)
subplot(2,1,1)
plot(GainFDFxyscan(:,1),GainFDFxyscan(:,2),'-oc','LineWidth',2);
% plot(wout/2/pi,(squeeze(mag)),'k','LineWidth',2);
hold on
plot(wout_FIR/2/pi,(squeeze(mag_FIR)),'k','LineWidth',2);
plot(wout_BJ/2/pi,(squeeze(mag_BJ)),'b','LineWidth',2);
plot(wout_FIR/2/pi,(squeeze(mag_FIR+1.96*sdmag_FIR)),'--k','LineWidth',1);
plot(wout_FIR/2/pi,(squeeze(mag_FIR-1.96*sdmag_FIR)),'--k','LineWidth',1);
% 
plot(wout_BJ/2/pi,(squeeze(mag_BJ+1.96*sdmag_BJ)),':b','LineWidth',2);
plot(wout_BJ/2/pi,(squeeze(mag_BJ-1.96*sdmag_BJ)),':b','LineWidth',2);

e = errorbar(GainFDFxyscan(:,1),GainFDFxyscan(:,2),err_g);
e.LineWidth = err_lw;
e.Color = err_c;
e.Marker = 'none';
% hold off
% legend('LES','Experimental','Harmonic','95% Confidence interval')
ax = gca;
ax.FontSize = 12;
ylabel ('$\mid{G}\mid$','FontSize',15,'Interpreter','latex')
axis([0 500 0 2])
grid on
set(gca,'xtick',[])
set(gca,'xticklabel',[])
xticks(0:100:500)
yticks(0:0.5:2)

% %--------Phase--------%
subplot(2,1,2)
phase0_FIR = phase_FIR(1);
phase0_BJ = phase_BJ(1);
h1 = plot(PhaseFDFxyscan(:,1),PhaseFDFxyscan(:,2).*-1,'-oc','LineWidth',2);
% h1 = plot(wout/2/pi,(squeeze(phase)-phase0)/180*pi,'k','LineWidth',2);
hold on
h2 = plot(wout_FIR/2/pi,(squeeze(phase_FIR)-phase0_FIR)/180*pi,'k','LineWidth',2);
h3 = plot(wout_BJ/2/pi,(squeeze(phase_BJ)-phase0_BJ)/180*pi,'b','LineWidth',2);


plot(wout_FIR/2/pi,(squeeze(phase_FIR+1.96*sdphase_FIR)-phase0_FIR)/180*pi,'--k','LineWidth',1);
plot(wout_FIR/2/pi,(squeeze(phase_FIR-1.96*sdphase_FIR)-phase0_FIR)/180*pi,'--k','LineWidth',1);

plot(wout_BJ/2/pi,(squeeze(phase_BJ+1.96*sdphase_BJ)-phase0_BJ)/180*pi,':r','LineWidth',2);
plot(wout_BJ/2/pi,(squeeze(phase_BJ-1.96*sdphase_BJ)-phase0_BJ)/180*pi,':r','LineWidth',2);

e = errorbar(PhaseFDFxyscan(:,1),-PhaseFDFxyscan(:,2),err_ph);
e.LineWidth = err_lw;
e.Color = err_c;
e.Marker = 'none';
hold off

legend([h1 h2 h3],'Experiment', 'FIR model','Box-Jenkins', 'Orientation','horizontal');
ax = gca;
ax.FontSize = 12;
ylabel ('$\angle{G}$','FontSize',15,'Interpreter','latex')
xlabel('Frequency [Hz]','FontSize',15,'Interpreter','latex')
axis([0 500 -16 4]) %%%%%%%%%%%%%%%%%% modif (-16)
grid on
xticks(0:100:500)
yticks(-16:4:4) %%%%%%%%%%%%%%%%%% modif (-16)
