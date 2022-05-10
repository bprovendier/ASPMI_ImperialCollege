%% 2.3. Adaptive Noise Cancellation %%

%% (a) Implementing ALE
clc; clear variables; close all;

M = [5:5:20]; delta = [1:25]; lr = 0.01;
eta = filter([1, 0, 0.5], [1], randn(1500, 1));
eta = eta(501:end);
x = sin([1:1000]*0.01*pi)';
s = x + eta;

MSPE = zeros(length(M), length(delta));
L = length(s); leg = {};

figure(1);
subplot(1, 2, 1); hold on; 
set(gca,'fontsize', 16); 
for i = 1:length(M)
    for j = 1:length(delta)
        [~,xhat,~] = ale_lms(s, lr, delta(j), M(i));
        MSPE(i, j) = (1/L)*(x-xhat)' * (x-xhat);
    end
    plot(delta, 10*log10(MSPE(i,:)));
    leg{i} = ['M = ' num2str(M(i))];
end
legend(leg); xlim([1, 25]);
xlabel('Delay ($\Delta$)', 'Interpreter', 'Latex');
ylabel('MSPE (dB)'); title('MSPE vs. Filter Order & Delay');

[w, xhat, error] = ale_lms(s, 0.01, 3, 5);
subplot(1, 2, 2); hold on;
set(gca,'fontsize', 16); 
plot([1:length(s)], s);
plot([1:length(xhat)], xhat);
plot([1:length(x)], x, 'LineWidth', 2); 
legend('Noisy signal', 'ALE', 'Clean signal')
xlabel('Sample (n)'); ylabel('Signal'); title('ALE Signal Denoising');
hold off;

%% (b) Testing out ALE system
clc; clear variables; close all;

lr = 0.01;
eta = filter([1, 0, 0.5], [1], randn(1500, 1));
eta = eta(501:end);
x = sin([1:1000]*0.01*pi)';
s = x + eta;
L = length(s);

M = 5; delta = 3;
MSPE = zeros(2, 2);

figure(1);
for i = 1:20
    [~,xhat,~] = ale_lms(s, lr, delta, i);
    MSPE(1, i) = (1/L)*(x-xhat)' * (x-xhat);
    [~,xhat,~] = ale_lms(s, lr, i, M);
    MSPE(2, i) = (1/L)*(x-xhat)' * (x-xhat);
end
subplot(1, 2, 1); hold on; set(gca,'fontsize', 16);
plot([1:20], 10*log10(MSPE(1,:)));  hold off;
xlabel('Filter Order (M)'); ylabel('MSPE (dB)');
title('MSPE vs. Filter Order ($\Delta=3$)', 'Interpreter', 'Latex'); 

subplot(1, 2, 2); hold on; set(gca,'fontsize', 16); 
plot([1:20], 10*log10(MSPE(2,:)),'r'); hold off;
xlabel('Delay ($\Delta$)', 'Interpreter', 'Latex'); 
ylabel('MSPE (dB)');
title('MSPE vs. Delay (M=5)', 'Interpreter', 'Latex'); 

%% (c) Testing out ANC system
clc; clear variables; close all;
lr = 0.01;
eta = filter([1, 0, 0.5], [1], randn(1500, 1));
eta = eta(501:end);
x = sin([1:1000]*0.01*pi)';
s = x + eta;
sec_noise = 0.5*eta - 0.1*delayseq(eta, 2); % delay of 2
L = length(s);

MSPE = zeros(2, 2);

[~,xhat] = ale_lms(s, lr, 3, 5);
[~,xhat_anc] = anc_lms(s, sec_noise, lr, 10);

figure(1); subplot(2, 1, 1);
hold on; set(gca,'fontsize', 18);
plot([1:length(xhat)], xhat);
plot([1:length(xhat_anc)], xhat_anc);
plot([1:length(x)], x, 'LineWidth', 2);
xlabel('Sample (n)'); ylabel('Signal'); 
title('ALE vs ANC For Signal Reconstruction');
legend('ALE', 'ANC', 'Clean');
hold off;

MSPE_ale = 10*log10(movmean((x-xhat).^2, 30));
MSPE_anc = 10*log10(movmean((x-xhat_anc).^2, 30));

subplot(2, 1, 2);
hold on; set(gca,'fontsize', 18);
plot([1:length(MSPE_ale)], MSPE_ale);
plot([1:length(MSPE_anc)], MSPE_anc);
xlabel('Sample (n)'); ylabel('MSPE (dB)'); 
title('MSPE: ALE vs ANC');
legend('ALE', 'ANC');
hold off;

%% (d)
clc; clear variables; close all;

load EEG_Data_Assignment2;
data = detrend(Cz);

t = [0:length(data)-1];
lr = [0.1, 0.01, 0.001];
M = [2, 5, 20];
stdevs = [2];
for i = stdevs
    mains = sin((2*pi*50/fs)*t) + (10^(-i))*randn(1, length(data));
    mains = mains';
    counter = 1;
    figure(i);
    sprintf('Variance: 1e-%d',i)
    for j = 1:length(lr)
        for k = 1:length(M)
            [w,xhat] = anc_lms(Cz, mains, lr(j), M(k));
            subplot(length(lr), length(M), counter)
            ylim([0, 55]); xlim([1, 4.8]); hold on; set(gca,'fontsize', 14);
            spectrogram(xhat, rectwin(3540), round(0.3*(3540)), 16382, fs, 'yaxis');
            title(['M = ',num2str(M(k)),' and \mu = ', num2str(lr(j))]);
            hold on;
            counter = counter + 1;
        end
    end
end

%% functions

function [w, xhat, error] = ale_lms(signal, lr, delta, M)

    w = zeros(M, length(signal)+1); 
    error = zeros(size(signal));
    xhat = signal(size(signal));
    
    for n = delta+M:length(signal)
        u = flip(signal(n-delta-M+1:n-delta));
        xhat(n) = dot(w(:, n), u);
        error(n) = signal(n)  - xhat(n);
        w(:, n+1) = w(:, n) + lr*error(n)*u;
    end
end

function [w, xhat] = anc_lms(signal, sec_noise, lr, M)

    w = zeros(M, length(signal)); 
    eta = zeros(size(signal));
    xhat = zeros(size(signal));
    u = delayseq(repmat(sec_noise, 1, M), [0:M-1])';
 
    for n = 1:length(signal)
        eta(n) = dot(w(:, n), u(:, n));
        xhat(n) = signal(n)  - eta(n);
        w(:, n+1) = w(:, n) + lr*xhat(n)*u(:, n);
    end
end