%% --- 4. From LMS to Deep Learning --- %%
%% (1) 
clc; clear variables; close all;
load('time-series.mat');

y = y - mean(y); 
order = 4;
x = [zeros(order,1); y]; 
[pr,error,y_hat] = lms_learning(y, x, order, 1e-5, 0, 'standard',0,false);

figure(1); 
subplot(1,2,1); hold on; set(gca,'fontsize', 16);
plot(1:length(y), y);
plot(1:length(y_hat), y_hat);
xlabel('Samples (n)');
ylabel('Y'); title('AR(4) Prediction');
legend('True Signal', 'AR(4)'); hold off;
subplot(1,2,2); hold on; set(gca,'fontsize', 16);
plot(1:length(y), y);
plot(1:length(y_hat), y_hat);
xlim([800,1000]);
xlabel('Samples (n)');
ylabel('Y'); title('AR(4) Prediction (Segment)');
legend('True Signal', 'AR(4)'); hold off;

MSE = (error'*error)/1000;
Rp = pow2db(var(y_hat)/var(error));
fprintf(sprintf('MSE: %.4f, Rp: %.4f',MSE,Rp));

%% (2 & 3) 
clc; clear variables; close all;
load('time-series.mat');

y = y - mean(y); 
order = 4;
x = [zeros(order,1); y]; 
[pr,error,y_hat] = lms_learning(y, x, order, 2e-7, 0, 'nonlinear',0.1,false);

figure(1); 
subplot(1,2,1); hold on; set(gca,'fontsize', 16);
plot(1:length(y), y);
plot(1:length(y_hat), y_hat);
xlabel('Samples (n)');
ylabel('Y'); title('Scaled LMS-Tanh Prediction');
legend('True Signal', 'LMS + tanh'); hold on;
subplot(1,2,2); hold on; set(gca,'fontsize', 16);
plot(1:length(y), y);
plot(1:length(y_hat), y_hat);
xlim([800,1000]);
xlabel('Samples (n)');
ylabel('Y'); title(' Scaled LMS-Tanh Prediction (Segment)');
legend('True Signal', 'LMS + tanh'); hold on;

MSE = (error'*error)/1000;
Rp = pow2db(var(y_hat)/var(error));
fprintf(sprintf('MSE: %.4f, Rp: %.4f',MSE,Rp));

%% (4) adding bias
clc; clear variables; close all;
load('time-series.mat');

order = 4;
x = [zeros(order,1); y];
[params,error,y_hat] = lms_learning(y, x, order, 2e-7, 0, 'nonlinear',1e-1,false);

figure(1); 
subplot(1,2,1); hold on; set(gca,'fontsize', 16);
plot(1:length(y), y);
plot(1:length(y_hat), y_hat);
xlabel('Samples (n)');
ylabel('Y'); title('Prediction with bias');
legend('True Signal', 'LMS + tanh'); hold off;
subplot(1,2,2); hold on; set(gca,'fontsize', 16);
plot(1:length(y), y);
plot(1:length(y_hat), y_hat);
xlim([800,1000]);
xlabel('Samples (n)');
ylabel('Y'); title('Prediction with bias (Segment)');
legend('True Signal', 'LMS + tanh'); hold off;

MSE = (error'*error)/1000;
Rp = pow2db(var(y_hat)/var(error));
fprintf(sprintf('MSE: %.4f, Rp: %.4f',MSE,Rp));

%% (5) Introducted pretraining
clc; clear variables; close all;
load('time-series.mat');

order = 4;
x = [zeros(order,1); y];
[params,error,y_hat] = lms_learning(y, x, order, 2e-7, 0, 'nonlinear',1e-1,true);

figure(1); 
subplot(1,2,1); hold on; set(gca,'fontsize', 16);
plot(1:length(y), y);
plot(1:length(y_hat), y_hat);
xlabel('Samples (n)');
ylabel('Y'); title('Prediction with pretraining');
legend('True Signal', 'LMS + tanh'); hold off;
subplot(1,2,2); hold on; set(gca,'fontsize', 16);
plot(1:length(y), y);
plot(1:length(y_hat), y_hat);
xlim([800,1000]);
xlabel('Samples (n)');
ylabel('Y'); title('Prediction with pretraining (Segment)');
legend('True Signal', 'LMS + tanh'); hold off;

MSE = (error'*error)/1000;
Rp = pow2db(var(y_hat)/var(error));
fprintf(sprintf('MSE: %.4f, Rp: %.4f',MSE,Rp));

%% functions

function [params, err, y_hat] = lms_learning(y, x, order, lr, gamma, opt, lr_a, pt)
    N = length(y);
    err = zeros(N,1);
    y_hat = zeros(N,1);

    if pt
        [params, a] = pretrain(y,x,100,20,order,lr,gamma,lr_a,opt);
    else
         params = zeros(order+1,N); 
         a = 30;
    end
    for n = 1:N
        if strcmp(opt,'standard')
            y_hat(n) = params(:,n).'*[1;x(n+order-1:-1:n)];
            act_factor = 1;
        elseif strcmp(opt,'nonlinear')
            nonlinear = tanh(params(:,n).'*[1;x(n+order-1:-1:n)]);
            y_hat(n) = a*nonlinear;
            act_factor = a*(1-(nonlinear.^2));
        else
            error("Please enter a valid algorithm option");
        end
        err(n) = y(n) - y_hat(n);
        if n < N
            params(:, n+1) = (1-lr*gamma)*params(:, n)...
                + lr*err(n)*act_factor*[1;x(n+order-1:-1:n)];
            if strcmp(opt,'nonlinear')
                a = a + lr_a*err(n)*nonlinear;
            end
        end
    end
end

function [params, a] = pretrain(y,x,n_epochs,K,order,lr,gamma,lr_a,opt)
    N = length(y);
    params = zeros(order+1,N);
    weights = zeros(order+1,1);
    err = zeros(N,1);
    y_hat = zeros(N,1);
    a = 40;
    for e = 1:n_epochs
        for k = 1:K
            if strcmp(opt,'standard')
                y_hat(k) = weights.'*[1;x(k+order-1:-1:k)];
                act_factor = 1;
            elseif strcmp(opt,'nonlinear')
                nonlinear = tanh(weights.'*[1;x(k+order-1:-1:k)]);
                y_hat(k) = a*nonlinear;
                act_factor = a*(1-(nonlinear.^2));
            else
                error("Please enter a valid algorithm option");
            end
            err(k) = y(k) - y_hat(k);
            if k < K
                weights = (1-lr*gamma)*weights...
                    + lr*err(k)*act_factor*[1;x(k+order-1:-1:k)];
                if strcmp(opt,'nonlinear')
                    a = a + lr_a*err(k)*nonlinear;
                end
            end
        end
    end
    params(:,1:K) = repmat(weights, 1, K);
end