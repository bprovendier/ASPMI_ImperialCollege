close all
clear all
clc

N=20;
n = 0:N-1;
nFFT = 256;
n_iter = 500;
Pseudospectrum = zeros(n_iter,nFFT);
m=(N/2)-1;

figure(1);
subplot(1,2,1)
for i=1:n_iter
    
    noise = 0.05/sqrt(2)*(randn(size(n))+1j*randn(size(n)));
    x = exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.32*n)+ noise;

    [X,R] = corrmtx(x,m,'mod');
    [S,F] = pmusic(R,2,nFFT,1,'corr');
    
    Pseudospectrum(i,:)=S;
    plot(F,S,'c'); set(gca,'xlim',[0.25 0.40]);
    hold on
end

avg_val = mean(Pseudospectrum);
plot(F,avg_val,'Linewidth',1);
grid on; xlabel('Normalised frequency (\pi rad/sample)','FontSize',9); ylabel('Pseudospectrum','FontSize',11);
title('MUSIC pseudospectrum','FontSize',11)

std_val = std(Pseudospectrum);
figure(1);
subplot(1,2,2)
plot(F,std_val,'Linewidth',1);
grid on; xlabel('Normalised frequency (\pi rad/sample)','FontSize',9); ylabel('Std pseudospectrum','FontSize',11);
set(gca,'xlim',[0.25 0.40]);
title('Standard Deviation of PSD estimate','FontSize',11)