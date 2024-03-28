clear;
clc;
close all
load orig11k.mat
% The script will first play back the given x11k at fs
% and then play the processed audio signal x8k820 for comparison.

%% Upsampling 
L = 4; % Upsampling factor
% Zero's padding/inserting
x44k1 = zeros(1,L*length(x11k));
x44k1(1:L:length(x44k1)) = x11k;

% LPF for upsampling:
%   Cutoff freq: wc = pi / L;
%   Gain L
P= 65;
h_LPF_up = L * FirFilterDesign(L,P);

% interpolate expanded signal
xi44k1 = filter(h_LPF_up,1,x44k1);

%% Downsampling:
M = 5;
x8k820 = xi44k1(1:M:length(xi44k1));

%% Playback
% play original sound
soundsc(x11k, fs);
pause(3.5)
% play processed sound
soundsc(x8k820, fs*4/5);
function hlpf = FirFilterDesign(L,P)
alpha = 1/L;
hlpf = alpha * sinc(alpha * (-(P-1)/2:(P-1)/2)).* hamming(P)';
end