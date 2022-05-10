%% --- 3. Adaptive signal processing --- %%

%% 3.2. Adaptive AR Model Based Time-Frequency Estimation %%

%% (a) 
clc; clear variables; close all;

n_samples = 1500;
f = [100*ones(500,1); ...
    100*ones(500,1) + ([501:1000]' - 500)/2; ...
    100*ones(500,1) + (([1001:1500]' -  1000)/25).^2];
phi = cumsum(f);
fs = 1500;
variance = 0.05;
eta = wgn(n_samples,1,pow2db(variance),'complex');
y = exp(1i*2*pi*phi/fs) + eta;
figure(1); subplot(1,2,1); hold on; set(gca,'fontsize', 16);
plot([1:n_samples], f);
xlabel('Time (samples)'); ylabel('Frequency (Hz)');
title('Frequency modulated (FM) signal'); hold off;
order = [1, 5, 10, 20];
subplot(1,2,2); hold on; set(gca,'fontsize', 16);
H = zeros(4,1500); w = zeros(4,1500);
for i = 1:4
    a = aryule(y,order(i));
    [H(i,:),w(i,:)] = freqz(1,a,n_samples,fs);
    P = abs(H(i,:)).^2;
    plot(w(i,:), pow2db(P));
end
plot(mean(f)*ones(1,1500), linspace(-10,30,1500),'--')
title('AR Spectrum estimates')
xlabel('Frequency (Hz)'); 
ylabel('Power/frequency (dB/(rad/sample))');
leg = arrayfun(@(a)strcat('p=',num2str(a)),order,'uni',0);
leg{end+1} = 'FMavg';
legend(leg); 
hold off;
%% (b)
clc; clear variables; close all;

n_samples = 1500;
f = [100*ones(500,1); ...
    100*ones(500,1) + ([501:1000]' - 500)/2; ...
    100*ones(500,1) + (([1001:1500]' -  1000)/25).^2];
phi = cumsum(f);
fs = 1500;
variance = 0.05;
eta = wgn(n_samples,1,pow2db(variance),'complex');
y = exp(1i*2*pi*phi/fs) + eta;

x = delayseq(y,1);
lrs = [1e-1, 1e-2,1e-3];
figure(1);
for idx = 1:3
    [a,~,~] = clms(y, x, 0, lrs(idx), 0);
    H = zeros(n_samples,n_samples);
    for n = 1:n_samples
        [h, w] = freqz(1 , [1; -conj(a(n))], n_samples,fs); 
        H(:, n) = abs(h).^2; 
    end

    medianH = 50*median(median(H));
    H(H > medianH) = medianH;
    subplot(3,1,idx); hold on; set(gca,'fontsize', 16);
    mesh([1:n_samples], w, H);
    view(2);
    xlabel('Time (samples)');
    ylabel('Frequency (Hz)');
    ylim([0,500]);
    title(strcat('CLMS Spectral Estimate ($\mu=',num2str(lrs(idx)),'$)'), 'Interpreter', 'Latex');
end

%% functions

function [h, error, y_hat] = clms(y, x, M, lr, leak)
    h = zeros(M+1, length(x),'like',1i); 
    error = zeros(size(x),'like',1i);
    y_hat = zeros(size(x),'like',1i);
    x_pad = [zeros(M,1); x];
    for n = 1:length(x)
        y_hat(n) = h(:,n)'*x_pad(n+M:-1:n); 
        error(n) = y(n) - y_hat(n);
        if n < length(x)
            h(:, n+1) = (1-lr*leak)*h(:, n) + lr*conj(error(n))*x_pad(n+M:-1:n);
        end
    end
end