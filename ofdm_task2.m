%%
% Practice2 Task2: OFDM channel estimation based on different pilot
% sequences

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
pilot_1 = ones([1, 17]) * 15;   
% 17 barker Pilot
pilot_2 = [-1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, -1,1,-1] * 15;
% 17 random Pilot
pilot_3 = rand([1, 17]) * 15;
pilot_group = {pilot_1, pilot_2, pilot_3};

interval = floor(256 / (length(pilot_1)-1)) + 1;
pilot_len = length(pilot_1);


%% OFDM

MMSE = zeros([3, 4]);  % store MMSE of 0dB, 10dB, 20dB, 30dB

for pilot_flag = 1:3
    
% construct Twiddle factor of [S] for linear model
pilot = pilot_group{pilot_flag};
N = 256 + pilot_len ;
const = exp(-1i*2*pi/N);
S_twiddle = zeros([pilot_len, 10]);
for k=0: pilot_len-1
    for l =0: 9
        
        S_twiddle(k+1, l+1) = const^(k * interval * l); % Twiddle factor
        
    end
end
S = S_twiddle .* repmat( pilot' , 1 , 10 );

for SNR = 0:1:30
    
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
MMSE(pilot_flag, SNR+1) = MSE_sum / times;

end
end

%% Plot
hold on;
plot(MMSE(1,1:end), 'k-');
plot(MMSE(2,1:end), 'k--');
plot(MMSE(3,1:end), 'k-.');
title('\bfComparison of different pilot sequences')
ylabel('MMSE')
xlabel('SNR (dB)')
legend('Optimal', 'Barker', 'Random');
grid on;
