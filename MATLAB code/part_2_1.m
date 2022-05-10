%% 2.1. The Least Mean Square (LMS) Algorithm %%
clc; clear all; close all;

%% (b)
clc; clear all; close all;
std_eta = 0.5;
data = filter([1], [1, -0.1, -0.8], std_eta*randn(1500, 100));
data = data(501:end, :);

% 1 realisation
[params1, error1] = lms_estimator(data(:, 1), 2, 0.05);
[params2, error2] = lms_estimator(data(:, 1), 2, 0.01);
figure(1); subplot(1, 2, 1); hold on; set(gca,'fontsize', 14);
plot([1:length(error1)], 10*log10(error1.^2));
plot([1:length(error2)], 10*log10(error2.^2));
xlabel('Sample'); ylabel('Squared Error (dB)');
title('Learning Curve (1 realisation)');
legend('$\mu=0.05$','$\mu=0.01$', 'Interpreter', 'Latex');
hold off;

% 100 realisations
error_tot = zeros(2, 100, 1000); lrs = [0.05, 0.01];
for i = 1:100
    for j = 1:2
        [a, error_tot(j, i, :)] = lms_estimator(data(:, i), 2, lrs(j)); 
    end
end
error_tot = squeeze(mean(10*log10(error_tot.^2), 2));
figure(1); subplot(1, 2, 2); hold on; set(gca,'fontsize', 14);
plot([1:length(error_tot(1, :))], error_tot(1, :));
plot([1:length(error_tot(2, :))], error_tot(2, :));
xlabel('Sample'); ylabel('Squared Error (dB)');
title('Learning Curve (100 realisations)');
legend('$\mu=0.05$','$\mu=0.01$', 'Interpreter', 'Latex');
hold off;

%% (c) 
clc; clear all; close all;

std_eta = 0.5;
data = filter([1], [1, -0.1, -0.8], std_eta*randn(200500, 100));
data = data(501:end, :); 
error_tot = zeros(2, 100, 200000); lrs = [0.05, 0.01];
for i = 1:100
    for j = 1:2
        [a, error_tot(j, i, :)] = lms_estimator(data(:, i), 2, lrs(j)); 
    end
end
error_tot = squeeze(mean(error_tot(:, 100, 1000:end).^2, [3, 2]));
M1 = (error_tot(1) - std_eta.^2)/std_eta.^2;
M2 = (error_tot(2) - std_eta.^2)/std_eta.^2;

%% (d) 
clc; clear all; close all;

a1 = 0.1; a2 = 0.8; std_eta = 0.5;
data = filter([1], [1, -0.1, -0.8], std_eta*randn(2000, 100));
data = data(501:end, :); 

params_tot = zeros(2, 2, 100, 1500); lrs = [0.05, 0.01];
for i = 1:100
    for j = 1:2
        [params_tot(j, :, i, :), e] = lms_estimator(data(:, i), 2, lrs(j)); 
    end
end
params_tot = squeeze(mean(params_tot, 3));
params_tot1 = squeeze(params_tot(1,:,:,:));
params_tot2 = squeeze(params_tot(2,:,:,:));
figure(1); subplot(1, 2, 1); hold on; set(gca,'fontsize', 18); ylim([0, 1]); xlim([0, 1000]);
h(1)=plot([1:length(params_tot1(1, :))], params_tot1(1, :));
h(2)=plot([1:length(params_tot1(2, :))], params_tot1(2, :));
h(3)=plot([1:length(params_tot1(1, :))], a1*ones(1,length(params_tot1(1, :))), '--');
h(4)=plot([1:length(params_tot1(2, :))], a2*ones(1,length(params_tot1(2, :))), 'k--');
xlabel('Time Step'); ylabel('Magnitude');
title('LMS Coefficient Estimation ($\mu=0.05$)', 'Interpreter', 'Latex');
legend(h,'$a_{1}$', '$a_{2}$', '$\hat{a}_{1}$','$\hat{a}_{2}$', 'Interpreter', 'Latex');
hold off;

subplot(1, 2, 2); hold on; set(gca,'fontsize', 18); ylim([0, 1]); xlim([0, 1000]);
h(1)=plot([1:length(params_tot2(1, :))], params_tot2(1, :));
h(2)=plot([1:length(params_tot2(2, :))], params_tot2(2, :));
h(3)=plot([1:length(params_tot2(1, :))], a1*ones(1,length(params_tot2(1, :))), '--');
h(4)=plot([1:length(params_tot2(2, :))], a2*ones(1,length(params_tot2(2, :))), 'k--');
xlabel('Time Step'); ylabel('Magnitude');
title('LMS Coefficient Estimation ($\mu=0.01$)', 'Interpreter', 'Latex');
legend(h,'$a_{1}$', '$a_{2}$', '$\hat{a}_{1}$','$\hat{a}_{2}$', 'Interpreter', 'Latex');
hold off;


%% (f) Leaky LMS
clc; clear all; close all;

a1 = 0.1; a2 = 0.8; std_eta = 0.5;
data = filter([1], [1, -0.1, -0.8], std_eta*randn(2000, 100));
data = data(501:end, :); 
lrs = [0.05, 0.01]; gammas = [0.1, 0.5, 0.9];
params_tot = zeros(2, 1500, 100); idx = [1 2; 3 4; 5 6];
figure(1);
for j = 1:2
    for k = 1:3
        for i = 1:100
            [params_tot(:, :, i), e] = leaky_lms(data(:, i), 2, lrs(j), gamma(k)); 
        end
        params = squeeze(mean(params_tot, 3));
        subplot(3, 2, idx(k, j)); hold on;
        set(gca,'fontsize', 18); ylim([0, 1]);
        h(1)=plot([1:length(params(1, :))], params(1,:),'LineWidth', 1); 
        h(2)=plot([1:length(params(2, :))], params(2, :),'LineWidth', 1);
        h(3)=plot([1:length(params(1, :))], a1*ones(1,length(params(1, :))), '--');
        h(4)=plot([1:length(params(2, :))], a2*ones(1,length(params(2, :))), '--');
        title(sprintf('($\\mu=%.2f, \\gamma=%.2f$)',...
            lrs(j), gammas(k)), 'Interpreter', 'Latex');
        legend(h,'$\hat{a}_{1}$','$\hat{a}_{2}$', '$a_{1}$', '$a_{2}$', 'Interpreter', 'Latex');
        xlabel('Time Step'); ylabel('Estimate');
        hold off
    end
end

%% functions

function [params, error] = lms_estimator(data, order, lr)

params = zeros(order, length(data));
error = zeros(size(data)) ;

for i = order+1:length(data)
    current_error = data(i) - dot(flip(data(i-order:i-1)), params(:, i-1));
    error(i) = current_error;
    params(:, i) = params(:, i-1) + lr*(current_error)*flip(data(i-order:i-1));
end
end

function [params, error] = leaky_lms(data, order, lr, gamma)

params = zeros(order, length(data)); 
error = zeros(size(data)) ;

for i = order+1:length(data)
    current_error = data(i) - dot(flip(data(i-order:i-1)), params(:, i-1));
    error(i) = current_error;
    params(:, i) = (1-lr*gamma)*params(:, i-1) + lr*(current_error)*flip(data(i-order:i-1));
end
end