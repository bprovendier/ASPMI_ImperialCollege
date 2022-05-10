%% 2.2. Adaptive Step Sizes %%

%% (a) Adaptive step size algorithms
clc; clear all; close all;

wo = 0.9; std_eta = (0.5).^0.5;
input = std_eta*randn(1500, 1000);
output = filter([1, wo], [1], input);
output = output(501:end, :);
input = input(501:end, :);

rho = 0.001; alpha = 0.8; lr0 = [0.2, 0.1, 0];
gass = {'standard', 'standard', 'ben', 'af', 'mx'};

figure(1);
for k = 1:3
    params_tot = zeros(5, 1000, 1000);
    lrs = [0.01, 0.1, lr0(k), lr0(k), lr0(k)]; 
    for j = 1:5
        for i = 1:1000
            [a, params_tot(j, :, i), e] = ...
                lms_arma(output(:, i), input(:, i), 0, 1, lrs(j), gass{j}, rho, alpha, 0);
        end
        param_error = -(squeeze(mean(params_tot(j, :, :), 3) - wo));
        subplot(3, 1, k); hold on; set(gca,'fontsize', 18);
        plot([1:length(param_error)], param_error); ylim([0, 1]);
        hold off
    end
    subplot(3, 1, k); title(sprintf('GASS Weight Error Curves ($\\mu_{0}=%.1f$)',lr0(k)), 'Interpreter', 'Latex');
    legend('$\mu=0.01$','$\mu=0.1$', 'Benveniste', 'Ang \& Farhang', 'Matthews \& Xie', 'Interpreter', 'Latex'); 
    ylabel('Weight Error'); xlabel('Time Step');
end

%% (c) GNGD
clc; clear all; close all;

wo = 0.9; std_eta = (0.5).^0.5;
input = std_eta*randn(750, 10000);
output = filter([1, wo], [1], input);
output = output(501:end, :);
input = input(501:end, :);

rho_g = 0.005; rho_b = 0.002;
lr_g = 1; lr_b = 0.1;

figure(1);
params_tot_g = zeros(250, 10000); error_tot_g = zeros(250, 10000);
params_tot_b = zeros(250, 10000); error_tot_b = zeros(250, 10000);
for i = 1:10000
    if mod(i, 100) == 0
    end
    [~, params_tot_g(:, i), error_tot_g(:, i)] = ...
        gngd(output(:, i), input(:, i), 0, 1, lr_g, rho_g, 0);
    [~, params_tot_b(:, i), error_tot_b(:, i)] = ...
        lms_arma(output(:, i), input(:, i), 0, 1, lr_b,'ben', rho_b, 0, 0);
end
params_tot_g = -(squeeze(mean(params_tot_g, 2)) - wo);
params_tot_b = -(squeeze(mean(params_tot_b, 2)) - wo);
error_tot_g = squeeze(mean(mag2db(error_tot_g.^2), 2));
error_tot_b = squeeze(mean(mag2db(error_tot_b.^2), 2));

subplot(1,2,1); xlim([0 80]); set(gca,'fontsize', 16); hold on;
plot([1:length(params_tot_g)], params_tot_g);
plot([1:length(params_tot_b)], params_tot_b);
title('Weight Error Dynamics')
legend('GNGD', 'Benveniste');
ylabel('Weight Error'); xlabel('Time Step');
hold off;

subplot(1,2,2); set(gca,'fontsize', 16); ylim([-20, 0]); xlim([0, 50]);
hold on;
plot([1:length(error_tot_g)], error_tot_g);
plot([1:length(error_tot_b)], error_tot_b);
title('Squared Prediction Error Dynamics');
legend('GNGD', 'Benveniste');
ylabel('Squared Prediction Error (dB)'); xlabel('Time Step');
hold off;

%% functions

function [ar_params, ma_params, error] = lms_arma(output, input, p, q, lr, gass, rho, alpha, leak)

    params = zeros(p+q, length(output)); 
    phi = zeros(p+q, length(output));
    lrs = lr*ones(size(output));
    error = ones(size(output));
    
    for i = max([p,q])+1:length(output)-1
        aug_dat = [flip(output(i-p:i-1)), flip(input(i-q:i-1))]';
        error(i) = output(i) - dot(aug_dat, params(:, i));
        params(:, i+1) = (1-leak*lr)*params(:, i) + lrs(i)*(error(i))*aug_dat;
        if strcmp(gass, 'af')
            if i > max([p,q])+1
                phi(:, i) = alpha*phi(:, i-1) + error(i-1)*prev_aug_dat;
            end
            lrs(i+1) = lrs(i) + rho*error(i)*dot(aug_dat, phi(:, i));
        elseif strcmp(gass, 'ben')
            if i > max([p,q])+1
                phi(:, i) = (eye(length(prev_aug_dat))-lrs(i-1)*(prev_aug_dat(:)*prev_aug_dat(:).'))...
                    *phi(:, i-1) + error(i-1)*prev_aug_dat;
            end
            lrs(i+1) = lrs(i) + rho*error(i)*dot(aug_dat, phi(:, i));
        elseif strcmp(gass, 'mx')
            if i > max([p,q])+1
                phi(:, i) = error(i-1)*prev_aug_dat;
            end
            lrs(i+1) = lrs(i) + rho*error(i)*dot(aug_dat, phi(:, i));
        elseif strcmp(gass, 'standard')
        end
        prev_aug_dat = aug_dat;
    end
    ar_params = params(1:p, :); ma_params = params(p+1:q, :);
end

function [ar_params, ma_params, error] = gngd(output, input, p, q, lr, rho, leak)

    params = zeros(p+q, length(output));
    error = ones(size(output));
    reg = ones(size(output))/lr;
    
    for i = max([p,q])+1:length(output)-1
        aug_dat = [flip(output(i-p:i-1)), flip(input(i-q:i-1))]';
        error(i) = output(i) - dot(aug_dat, params(:, i));
        lr_now = lr/(reg(i) + dot(aug_dat, aug_dat));
        params(:, i+1) = (1-leak*lr)*params(:, i) + lr_now*(error(i))*aug_dat;
        if i >  max([p,q])+1
            num = rho*lr*error(i)*error(i-1)*dot(old_dat, aug_dat);
            den = ( reg(i-1) + dot(old_dat, old_dat) ).^2;
            reg(i+1) = reg(i) - num/den;
        end
        old_dat = aug_dat;
        
    end
    ar_params = params(1:p, :); ma_params = params(p+1:q, :);
end