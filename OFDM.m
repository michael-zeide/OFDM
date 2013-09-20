clear all;
close all;
clc;

%% Parameters

PPM = 20;
R = 1e5; % [bits/sec]
duration = 0.128; % [sec]
DataL = R*duration/2;  % Data length in symbols
NFFT = 128;

M = 4; % 4 QAM
snr = 20; % [dB]

Nsym = 4;           % Filter order in symbol durations
beta = 0.5;         % Roll-off factor
sampsPerSym = 10;    % Upsampling factor
L = sampsPerSym*Nsym + 1; % Raised cosine filter order
Fs = R * sampsPerSym;   % Sampling frequency
f1 = Fs/5;
cyclic_prefix_signal = [];

%% Raised cosine filter design

shape = 'Raised Cosine';
% Specifications of the raised cosine filter with given order in symbols
rcosSpec = fdesign.pulseshaping(sampsPerSym, shape, 'Nsym,beta', Nsym, beta);
rcosFlt = design(rcosSpec);
rcosFlt.Numerator = rcosFlt.Numerator / max(rcosFlt.Numerator);


%% Tranceiver

bernoulli_binary_generator = randint(R*duration,1);

% Transmiter

bernoulli_two_samples = reshape(bernoulli_binary_generator, length(bernoulli_binary_generator)/2, 2);
dec = bi2de(bernoulli_two_samples,'left-msb'); % Bit to integer
modulated_data = qammod(dec,M); % 4 QAM

% Serial to Parallel
for id = 1:length(modulated_data)/NFFT
    ifft_signal = ifft(modulated_data((id-1)*NFFT+1:id*NFFT)); % ifft
    % adding cyclic prefix
    cyclic_prefix = zeros(NFFT + L, 1);
    cyclic_prefix(1:L) = ifft_signal(end-L+1:end);
    cyclic_prefix(L+1:end) = ifft_signal;
    cyclic_prefix_signal = [cyclic_prefix_signal; cyclic_prefix];
end;  
% D/A
signal_complex = filter(rcosFlt, upsample([cyclic_prefix_signal; zeros(Nsym/2,1)], sampsPerSym));
fltDelay = Nsym / (2*R); % Filter group delay, since raised cosine filter is linear phase and symmetric.
signal_complex = signal_complex(fltDelay*Fs+1:end); % Correct for propagation delay by removing filter transients
    
t = (0: length(signal_complex) - 1) / Fs;
f = linspace(-Fs/2,Fs/2,length(signal_complex));
I = real(signal_complex);
Q = imag(signal_complex);
FI = fftshift(fft(I));
FQ = fftshift(fft(Q));
signal = I'.*cos(2*pi*f1*t) - Q'.*sin(2*pi*f1*t);
    
% Channel    
noised_signal = awgn(signal,snr,'measured'); % Adding white gaussian noise
% noised_signal = signal;
    
% Reciever
    
% A/D
noised_recieved_signal = noised_signal.*cos(2*pi*f1*t) - 1i*noised_signal.*sin(2*pi*f1*t);

filtered_reconstructed_digital_signal = filter(rcosFlt, [noised_recieved_signal.'; zeros(fltDelay*Fs,1)]); 
reconstructed_digital_signal = downsample(filtered_reconstructed_digital_signal, sampsPerSym);
reconstructed_digital_signal = reconstructed_digital_signal(Nsym/2+1:end); % Correct for propagation delay by removing filter transients

% Serial to Parallel
N = L + NFFT;
fft_signal = zeros(length(modulated_data),1);
for id = 1:length(modulated_data)/NFFT
    tmp = reconstructed_digital_signal(N*(id-1) + 1 : N*id);
    removing_cyclic_prefix = tmp(end - NFFT + 1:end); % removing cyclic prefix
    fft_signal_SP = fft(removing_cyclic_prefix); % fft
    fft_signal(NFFT*(id-1) + 1 : id*NFFT) = fft_signal_SP;
end;


demodulated_data = qamdemod(fft_signal, M); % Demodulation
demodulated_binary = de2bi(demodulated_data,'left-msb'); % Integer to bit
recieved_signal = demodulated_binary(:);
    
I = real(fft_signal);
Q = imag(fft_signal);
scatter(I,Q);
hold on;
scatter([1,-1,-1,1],[1,1,-1,-1]);
title(sprintf('Reciever''s constelation, SNR = %d[dB]', snr));

BER = sum(xor(bernoulli_binary_generator,recieved_signal))/length(bernoulli_binary_generator);