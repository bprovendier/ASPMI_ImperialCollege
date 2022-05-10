close all
clear all
clc

%% General parameters
L=5000;
K=2000;
M=500; 
fs=80;
f = [0.9, 1.2, 1.5];
psd = zeros(M,2*L-1);
t_axis = (0:K-1)/fs;
freqAxis = (-(2*L-1)/2:(2*L-1)/2-1)*(fs/(2*L-1));

signal=[0.4*sin(2*pi*t_axis*f(1)) + 0.6*sin(2*pi*t_axis*f(2)) + 0.5*sin(2*pi*t_axis*f(3)) zeros(1, L-length(t_axis))];

figure(1);
subplot(1,2,1)
for i=1:M
    corrupted_signal = signal + randn(1,L);
    [acf,~] = xcorr(corrupted_signal, 'biased');
    psd(i,:) =real(fftshift(fft(ifftshift(acf))));
    plot(freqAxis,psd(i,:),'c','Linewidth',1);
    hold on
end
xlim([0 2])
avg_psd = mean(psd);
plot(freqAxis,avg_psd,'b','Linewidth',1);
xlabel('Frequency (Hz)','FontSize',11); ylabel('PSD','FontSize',11)
title({'PSD estimates','(realisations and mean)'},'FontSize',11)
std_psd = std(psd);

subplot(1,2,2)
plot(freqAxis,std_psd,'r','Linewidth',1);
xlabel('Frequency (Hz)','FontSize',11); ylabel('Standard deviation','FontSize',11)
title({'Standard deviation','of the PSD estimate'},'FontSize',11)
xlim([0 2])

figure(2);
subplot(1,2,1)
for i=1:M
   psd_dB(i,:) = pow2db(psd(i,:));
   plot(freqAxis,psd_dB(i,:),'c','Linewidth',1); 
   hold on 
end
xlim([0 2])
avg_psd_dB = mean(psd_dB);
plot(freqAxis,avg_psd_dB,'b','Linewidth',1);
xlabel('Frequency (Hz)','FontSize',11); ylabel('PSD (dB)','FontSize',11)
title({'PSD estimates in dB','(realisations and mean)'},'FontSize',11)

std_psd_dB = std(psd_dB);
subplot(1,2,2)
plot(freqAxis,std_psd_dB,'r','Linewidth',1);
xlim([0 2])
xlabel('Frequency (Hz)','FontSize',11); ylabel('Standard deviation (dB)','FontSize',11)
title({'Standard deviation','of the PSD estimate (dB)'},'FontSize',11)