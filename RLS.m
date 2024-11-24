

% Recursive Least Squares (RLS) Based on Provided Equations
clear; clc;

% Parameters
N = 120;                % Number of iterations (warming-up period)
lambda = 0.98;          % Forgetting factor
n_theta = 2;            % Number of parameters
P = eye(n_theta);       % Initial covariance matrix (identity matrix)
Theta = zeros(n_theta, 1); % Initial parameter estimates

% Synthetic Data Generation
x = randn(N, 1);        % Input data
y = zeros(N, 1);        % Output data
theta_true = [1.2; -0.8]; % True parameters of the system

for t = 2:N
    y(t) = theta_true(1) * x(t) + theta_true(2) * x(t-1) + 0.1 * randn; % Generate output
end

% Recursive Least Squares Implementation
Theta_history = zeros(N, n_theta); % Store parameter estimates
error_history = zeros(N, 1);       % Store errors

for k = 2:N
    % Create regression vector Z_k
    Z_k = [x(k); x(k-1)];
    
    % Update P_k+1
    K_k = (P * Z_k) / (lambda + Z_k' * P * Z_k); % Kalman gain
    P = (P / lambda) - K_k * Z_k' * P / lambda;
    
    % Update Theta_k+1
    prediction_error = y(k) - Z_k' * Theta;
    Theta = Theta + P * Z_k * prediction_error;

    % Store results
    Theta_history(k, :) = Theta';
    error_history(k) = prediction_error;
end

% Plot Results
figure;
subplot(2, 1, 1);
plot(Theta_history);
hold on;
plot(1:N, repmat(theta_true', N, 1), '--');
xlabel('Time Step');
ylabel('Parameter Estimates');
legend('\Theta_1', '\Theta_2', 'True \Theta_1', 'True \Theta_2');
title('RLS Parameter Estimates');

subplot(2, 1, 2);
plot(error_history);
xlabel('Time Step');
ylabel('Prediction Error');
title('Prediction Error Over Time');
