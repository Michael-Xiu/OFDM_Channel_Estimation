%%
% Practice2 Task1: OFDM channel estimation based on pilot sequence in
% different SNR situations

% Course: Random Signal Analysis, Professor Dai
% Author: Xiu Shengjie, 16308125
% School of Electronics and Information Technology, SYSU
% Date: 6/13/2019


%% Field
clear
close all


h=[1:-0.1:0.1];
times = 1000;

% 17 Unity Pilot
pilot = ones([1, 17]) * 15;   

interval = floor(256 / (length(pilot)-1)) + 1;
pilot_len = length(pilot);

% construct Twiddle factor of [S] for linear model
N = 256 + pilot_len ;
const = exp(-1i*2*pi/N);
S_twiddle = zeros([pilot_len, 10]);
for k=0: pilot_len-1
    for l =0: 9
        
        S_twiddle(k+1, l+1) = const^(k * interval * l); % Twiddle factor
        
    end
end
S = S_twiddle .* repmat( pilot' , 1 , 10 );


%% OFDM

MMSE = zeros([1, 4]);  % store MMSE of 0dB, 10dB, 20dB, 30dB

for SNR = 0:10:30
    
MSE_sum = 0;

for time = 1: times

%% Input: Gaussian
x_fft =wgn(1, 256, 0);

%% Pilot Insert
x_fft_pilot = zeros([1, 256+pilot_len]);
count_pilot = 1;
count_data = 1;
for k = 1:length(x_fft_pilot)
    if mod(k - 1, interval) == 0 && count_pilot <=pilot_len;
        x_fft_pilot(k) = pilot(count_pilot);
        count_pilot =count_pilot + 1;
    else
        x_fft_pilot(k) = x_fft(count_data);
        count_data =count_data + 1;
    end
end

%% IFFT
x_ifft_pilot=ifft(x_fft_pilot);

%% Adding Cyclic Extension
x_cp=zeros([1, length(x_ifft_pilot) + 10]);
x_cp(1:10) = x_ifft_pilot(end-9:end);
x_cp(11:end) = x_ifft_pilot;

%% Channel
y_cp = conv(x_cp, h);
y_cp = awgn(y_cp, SNR);

%% Removing Cyclic Extension and its tail
y_ifft = y_cp(11:end-9);

%% FFT
y_fft_pilot = fft(y_ifft); % retain the same FFT length as original x

%% Pilot Extract
pilot_Rx =  y_fft_pilot(1: interval: end);

%% Linear Model MVU Estimator
h_hat = (S' * S)^-1 * S' * pilot_Rx.';
h_hat = abs(h_hat);

%% Error: MSE
MSE = sum( (h-h_hat.').^2 ./ h.^2 );
MSE_sum = MSE_sum+MSE;

end
MMSE(SNR/10+1) = MSE_sum / times;

end


%% Plot
b=diag(MMSE);
c=bar(b,'stacked');
legend('0dB','10dB', '20dB', '30dB')
title('OFDM Channel Estimation Based on Linear Model')
ylabel('MMSE')
