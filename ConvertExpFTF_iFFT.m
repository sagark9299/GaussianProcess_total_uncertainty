clear
clc

% Load FTF data
% experiment = load('FTF_A.mat');
% load FTF_gain_10_new_arc.mat
% load FTF_phase_10_new_arc.mat
load experiment_app.mat
% experiment.gain = GainFDFxyscan(:,2);
% experiment.freq = GainFDFxyscan(:,1);
% experiment.phase = PhaseFDFxyscan(:,2);
N = 250;
experiment_app.phase = -experiment_app.phase;
% experiment.phase = -experiment.phase;
% Append the FTF data for conversion
% 1-interpolate gain/phase (200-400 Hz)
% Freq_interp = 200:10:400;
% gain_app = interp1(experiment.Freq(21:end),experiment.gain(21:end),Freq_interp,'spline');
% phase_app = interp1(experiment.Freq(21:end),experiment.phase(21:end),Freq_interp,'spline');
% 2-Append Freq/gain/phase (0-400 Hz)
% Freq_list = [experiment.Freq(1:20);Freq_interp'];
% gain_list = [experiment.gain(1:20);gain_app'];
% phase_list = [experiment.phase(1:20);phase_app'];
% 3-Append Freq/gain/phase (400-1250 Hz)
Freq_interp = (0:10:480)';
gain_list = experiment_app.gain;
phase_list = -experiment_app.phase;
freq_list = experiment_app.freq;
gain_half = [gain_list;zeros(N/2-size(gain_list,1)+1,1)];  % For gain, add zero
phase_coeff = polyfit(Freq_interp,phase_list,2);            % For phase, linear extrapolation
phase_extrap = polyval(phase_coeff,490:10:1250); % 410:10:1250
% phase_extrap = interp1(Freq_interp,phase_app,410:10:1250,'linear','extrap'); % For phase, linear extrapolation
phase_half = [phase_list;phase_extrap'];
% Generate conjugate
FTF_half = gain_half.*exp(1i*phase_half);
FTF_conj = conj(FTF_half);
FTF_full = [FTF_half;FTF_conj(end-1:-1:2)];

% Convert back to FIR model
complex_vector = zeros(N,N);
for index_i = 1:N
    for index_j = 1:N
        complex_vector(index_i,index_j) = exp(1i*2*pi/N*(index_i-1)*(index_j-1));
    end
end
FIR_coeff = real(1/N*FTF_full.'*complex_vector);

% Post-processing
% Impulse response
time = 0:4.1e-4:(N-1)*4.1e-4;
figure(1)
% stem(time(1:30),FIR_coeff(1:30),'k','filled','LineWidth',1.2,'MarkerSize',6);
plot(time(1:30),FIR_coeff(1:30),'k','LineWidth',1.2);
% title('Experimental FIR')
xlabel('Time (s)')
ylabel('Normalized FIR coefficient (-)')
ax = gca;
ax.FontSize = 12;


% Frequency response
Freq_evaluate = 0:10:480;
[gain, phase] = FTF_construct(FIR_coeff, 4e-4, Freq_evaluate');
figure(2)
subplot(2,1,1)
plot(Freq_evaluate,gain,'o','MarkerSize',8)
subplot(2,1,2)
plot(Freq_evaluate,phase,'o','MarkerSize',8)
axis([0 500 -16 4])

% save 'FIR_A.mat' FIR_coeff