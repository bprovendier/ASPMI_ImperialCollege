clear; close all; clc;
fSample = 1;
nSamples = 1e4;
nTransients = 5e2;
coefAr = [2.76 -3.81 2.65 -0.92];
variance = 1;
orderAr = 2: 14;
nOrders = length(orderAr);
arModel = arima('AR', coefAr, 'Variance', variance, 'Constant', 0);
arSignal = simulate(arModel, nSamples);
arSignal = arSignal(nTransients + 1: end);
nSamples = length(arSignal);
[h, f] = freqz(1, [1 -coefAr], nSamples, fSample);
psd = abs(h) .^ 2;
varEst = zeros(nOrders, 1);
psdAr = cell(nOrders, 1);
for iOrder = 1: nOrders
    [coefArEst, varEst(iOrder)] = aryule(arSignal, orderAr(iOrder));
    hAr = freqz(sqrt(varEst(iOrder)), coefArEst, nSamples);
    psdAr{iOrder} = abs(hAr) .^ 2;
end

figure;
% PSD: ground truth vs AR estimation
subplot(2, 1, 1);
plot(f, pow2db(psd), 'LineWidth', 2);
hold on;
plot(f, pow2db(psdAr{1}), 'LineWidth', 2);
hold on;
plot(f, pow2db(psdAr{7}), 'LineWidth', 2);
hold on;
plot(f, pow2db(psdAr{13}), 'LineWidth', 2);
grid on;
legend('Truth', 'AR (2)', 'AR (8)', 'AR (14)');
title(sprintf('PSD estimate by AR model: signal length = %d', nSamples));
xlabel('Normalised frequency (\pi rad/sample)');
ylabel('PSD (dB)');
% variance (noise power)
subplot(2, 1, 2);
plot(orderAr, pow2db(varEst), 'LineWidth', 2);
grid on;
legend('Error');
title('Prediction error against AR order');
xlabel('AR order');
ylabel('Noise power (dB)');