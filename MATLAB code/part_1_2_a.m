clear all
close all
clc

%% Load the sunspot data
load sunspot.dat
fs=1; %Sampling frequency - One value per year

time = sunspot(:,1);
sunspot_series = sunspot(:,2);
L=length(sunspot_series);

%% Adding small a DC values to avoid zero values in log
sunspot_series=sunspot_series+eps;%distance from 1.0 to the next largest double-precision numbe

[psd_original f] = periodogram(sunspot_series,hann(L),L,fs,'one-sided');

figure(4);
subplot(1,3,1)
plot(f,psd_original,'Linewidth',1)
ylabel('PSD','FontSize',11);
xlabel('Cycles/Year','FontSize',11);
title('PSD est. of pre-processed data')
xlim([0 fs/2])

%% Removing mean and trend of series

%Remove mean only
sunspot_zero_mean_env = sunspot_series - mean(sunspot_series);

%Removing linear trend as well
sunspot_detrend = detrend(sunspot_zero_mean_env,1);

fft_detrend = fftshift(fft(sunspot_detrend,L));

[psd_detrend,f] = periodogram(sunspot_detrend,hann(L),L,fs,'one-sided');

figure(4);
subplot(1,3,2)
plot(f,psd_detrend,'Linewidth',1)
ylabel('PSD','FontSize',11);
xlabel('Cycles/Year','FontSize',11);
title('PSD est. of centered&detrended data')
xlim([0 fs/2])

%% Considering logarithm of data and removing mean
sunspot_log = log(sunspot_series)-mean(log(sunspot_series));

%Considering modified periodogram method (periodogram-based technique)
[psd_estimate f] = periodogram(sunspot_log,hann(L),L,fs,'one-sided');

figure(4);
subplot(1,3,3)
plot(f,psd_estimate,'Linewidth',1);
ylabel('PSD','FontSize',11);
xlabel('Cycles/Year','FontSize',11);
title('PSD est. of centered logarithmic data')
xlim([0 fs/2])