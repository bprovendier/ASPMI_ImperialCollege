close all
clear all
clc
fs=1;
N=[20,40,110];

figure(1);
for i=1:3
    n = 0:N(i)-1;
    noise = 0.2/sqrt(2)*(randn(size(n))+1j*randn(size(n)));
    x = exp(1j*2*pi*0.3*n)+exp(1j*2*pi*0.32*n)+ noise;

    dF = fs/N(i); 
    dF_new=fs/512;
    K=dF_new/dF;
    [pxx,f] = periodogram(x,rectwin(N(i)),round(N(i)/K),fs);
    plot(f,pxx,'Linewidth',1)
    hold on
end
xlabel('Frequency (Hz)','FontSize',11)
ylabel('PSD','FontSize',11)
grid on
title('PDS estimates of complex exponentials','FontSize',11)
legend('N=20','N=40','N=110','FontSize',9)
xlim([0.2 0.4])
