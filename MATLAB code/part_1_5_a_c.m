%% Part a
close all; clear all; clc;

% Add data folder to path
addpath('C:\Users\bapti\OneDrive - Imperial College London\Documents\MATLAB\RRI_data');
load RRI-DATA.mat;

xRRI1 = detrend(normalize(xRRI1));
xRRI2 = detrend(normalize(xRRI2));
xRRI3 = detrend(normalize(xRRI3));
fsRRI1 =5;
RRI_data = {xRRI1; xRRI2; xRRI3};
leg = {'Standard', '(W = 50s)', '(W = 150s)'};
figure(1);
for j = 1:3 % Trials
    subplot(3, 1, j); hold on;
    L = length(RRI_data{j}); wl = [L, 50*fsRRI1, 150*fsRRI1];
    for i = 1:3 % Window sizes
        [P_mean, w] = pwelch(RRI_data{j}, wl(i), 0, L, 4);
        plot(w, 10*log10(P_mean));
    end
    set(gca,'fontsize', 14);
    xlabel('Frequency (Hz)');
    ylabel('Pow/freq (dB/Hz)'); 
    title(sprintf('RRI PSD Estimate (Trial %d)', j));
    legend(leg); hold off;
end

%% Part c
clc; clear all; close all;
% Add data folder to path
addpath('C:\Users\bapti\OneDrive - Imperial College London\Documents\MATLAB\RRI_data');
load RRI-DATA.mat;

xRRI1 = detrend(normalize(xRRI1));
xRRI2 = detrend(normalize(xRRI2));
xRRI3 = detrend(normalize(xRRI3));


order = [3, 7, 4]; RRI_data = {xRRI1; xRRI2; xRRI3};
figure(1); hold on;
for i = [1, 2, 3]
    [pxx, w] = pyulear(RRI_data{i}, order(i), 2048, 4); % power spectrum estimate given AR model
    plot(w, 10*log10(pxx));
end
set(gca,'fontsize', 14);
xlabel('Frequency (Hz)');
ylabel('Pow/freq (dB/Hz)'); 
title('AR PSD Estimate of RRI Data');
legend('Trial 1', 'Trial 2', 'Trial 3'); hold off;
