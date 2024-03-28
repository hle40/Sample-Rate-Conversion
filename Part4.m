clear;
clc;
close all
load orig11k.mat

N=2^12; % 2048
%% The original signal frequency domain
X11k = fft(x11k,N);

%% Upsampling 
L = 4; % Upsampling factor

% Zero's padding/inserting
x44k1 = zeros(1,L*length(x11k));
x44k1(1:L:length(x44k1)) = x11k;
% Compute fft of upsampled signal
X44k1 = fft(x44k1,N);

% LPF for upsampling:
%   Cutoff freq: wc = pi / L;
%   Gain L
P= 65;
h_LPF_up = L * FirFilterDesign(L,P);
H_LPF_up = fft(h_LPF_up,N);

% interpolate expanded signal
xi44k1 = filter(h_LPF_up,1,x44k1);
Xi44k1 = fft(xi44k1,N);

%% Downsampling:
M = 5;
x8k820 = xi44k1(1:M:length(xi44k1));
X8k820 = fft(x8k820,N);

% Filter for Downsampling:
%   Cutoff freq: wc = pi / M;
%   Gain 1
h_LPF_down = FirFilterDesign(M,P);
H_LPF_down = fft(h_LPF_down,N);

% anti alias filter prior to downsampling:
prior_down = filter(h_LPF_down,1,xi44k1);
X_prior_down = fft(prior_down,N);

% Downsampled the filtered signal:
x8k820_filtered = prior_down(1:M:length(xi44k1));
X8k820_filtered = fft(x8k820_filtered,N);

%% Plotting
% Original Signal
figure(1)
plotMagFreq(X11k,N)
ylabel('|X_{11k}(e^{j\omega})|')
title('DFT of Original Signal')

% Upsampled signal
figure(2)
plotMagFreq(X44k1,N)
ylabel('|X_{44k1}(e^{j\omega})|')
title('DFT of Upsampled Signal')

% LPF for upsampling
figure(3)
k=0:(N-1);
wk= 2*k/N;
plot(wk,angle(H_LPF_up),'LineWidth',1.5)
xlabel('\omega_\pi (DT freq normalized by \pi)')
ylabel('\angle H(e^{j\omega})')
xticks(0:0.1:2)
grid on
title('Phase of Frequency Response of LPF for Upsampling')

figure(4)
plotMagFreq(H_LPF_up,N)
ylabel('|H(e^{j\omega})|')
title('Frequency Repsonse of LPF for Upsampling')

% Plot interploted signal for upsampling:
figure(5)
plotMagFreq(Xi44k1,N)
ylabel('|X_{interpolated}(e^{j\omega})|')
title('DFT of filtered upsampled signal')
hold on
plotMagFreq(H_LPF_down,N)
legend('','LPF prior to Downsampling M=5')

% plot DFT of downsampled singal
figure(6)
plotMagFreq(X8k820,N)
ylabel('|X_{downsampled}(e^{j\omega})|')
title('DFT of downsampled signal (no prior filter)')
hold on

% plot filtered signal prior to downsampling
figure(7)
plotMagFreq(X_prior_down,N)
ylabel('|X_{prior down}(e^{j\omega})|')
title('DFT of signal LPFed in preparation for downsampling')
hold on

% DFT of filtered prior to downsampling
figure(8)
plotMagFreq(X8k820_filtered,N)
ylabel('|X_{8820}(e^{j\omega})|')
title('DFT of LPFed signal and then downsampled')

%saveAllFigure(8)
function plotMagFreq(X,N)
k=0:(N-1);
wk= 2*k/N;
plot(wk,abs(X),'LineWidth',1.5)
xlabel('\omega_\pi (DT freq normalized by \pi)')
xticks(0:0.1:2)
grid on
end

function stemPlot(n,x)
stem(n,x)
xlabel('n')
grid on
end

function hlpf = FirFilterDesign(L,P)
alpha = 1/L;
hlpf = alpha * sinc(alpha * (-(P-1)/2:(P-1)/2)).* hamming(P)';
end
