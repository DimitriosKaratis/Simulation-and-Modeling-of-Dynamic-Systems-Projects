%% SIMULATION AND MODELING OF DYNAMIC SYSTEMS
% Final Project – Exercise 2b
% Dimitrios Karatis – 10775

clear; clc; close all;
format long

global c d f gamma thetam theta_true newInput

% True system parameters
% The true function is: f(x, u, θ) = -x^3 + θ1 * tanh(x) + θ2 / (1 + x^2) + u
theta_true = [0.8, 1.2];

% Estimator gains
thetam = 41.2;
gamma = [39.9 41.7 43.3];

% Simulation settings
tStart = 0; tEnd = 200; tStep = 0.01;
tspan = tStart:tStep:tEnd;
N = length(tspan);

% Original input signal u = c * sin(d * t) + f + 0.1 * sin(0.2*t) + 0.5 * sin(0.3*t) parameters
c = 2.5; d = 1; f = 0;
newInput = false;  % Use original input first

% Storage for results
energy_error = zeros(1,5);
BIC = zeros(1,5);
k_vec = [2 3 3 3 3];  % Number of parameters per model

% === Model evaluation loop (original input) ===
for modelType = [1 2 3 4 5]
    initCond = [0; 0; 0.01; 0.01; 0.01];  % Initial conditions
    [t,x] = ode45(@(tpar,xpar) LyapSP(tpar,xpar,modelType), tspan, initCond);

    e1 = x(:,1) - x(:,2);                   % Tracking error
    E = trapz(t, e1.^2);                    % Total error energy
    energy_error(modelType) = E;

    MSE = E / N;                            % Mean square error
    k = k_vec(modelType);                   % Number of parameters
    BIC(modelType) = k*log(N) + N*log(MSE); % Bayesian Information Criterion
end

% --- Print Results ---
fprintf('\n=== BIC Analysis ===\n');
for i = 1:5
    fprintf('Model %d | E = %.5f | BIC = %.2f\n', i, energy_error(i), BIC(i));
end

BIC(3) = Inf;  % Exclude model 3 from selection
[~, bestModel] = min(BIC);
fprintf('\nBest model (excluding M3): Model %d\n', bestModel);


%% === Stability Test with New Input ===
newInput = true;  % Use new input signal
fprintf('\n=== Stability Test with New Input ===\n');

initCond = [0; 0; 0.01; 0.01; 0.01];
[t2,x2] = ode45(@(t,x) LyapSP(t,x,bestModel), tspan, initCond);

e2 = x2(:,1) - x2(:,2);
theta1_hat2 = x2(:,3); theta2_hat2 = x2(:,4); theta3_hat2 = x2(:,5);
E2 = trapz(t2, e2.^2);
fprintf('Error energy with new input: %.7f\n', E2);

% === Re-run with original input for comparison ===
newInput = false;
[t1,x1] = ode45(@(t,x) LyapSP(t,x,bestModel), tspan, initCond);
e1 = x1(:,1) - x1(:,2);
theta1_hat1 = x1(:,3); theta2_hat1 = x1(:,4); theta3_hat1 = x1(:,5);

% === Plot comparison ===
figure;

subplot(3,1,1);
plot(t1, x1(:,1), 'r', t1, x1(:,2), 'b--'); hold on;
plot(t2, x2(:,1), 'm', 'LineWidth', 2.5); 
plot(t2, x2(:,2), 'c--', 'LineWidth', 2.5);
legend('x_1 (orig)', 'x1_{hat} (orig)', 'x1 (new)', 'x1_{hat} (new)');
title('Comparison of x_1 and x1_{hat} under original and new input');
xlabel('Time (t)'); ylabel('State'); grid on;

subplot(3,1,2);
plot(t1, e1, 'k', 'LineWidth', 1); hold on;
plot(t2, e2, 'g--', 'LineWidth', 1);
legend('e_1 (orig input)', 'e_1 (new input)');
title('Tracking Error Comparison');
xlabel('Time (t)'); ylabel('e_1(t)'); grid on;

subplot(3,1,3);
plot(t1, theta1_hat1, 'r'); hold on;
plot(t1, theta2_hat1, 'b');
plot(t1, theta3_hat1, 'k');
plot(t2, theta1_hat2, 'm--', 'LineWidth', 1.5);
plot(t2, theta2_hat2, 'c--', 'LineWidth', 1.5);
plot(t2, theta3_hat2, 'g--', 'LineWidth', 1.5);
title('Parameter Estimates θi_{hat}(t)');
legend('θ1_{hat} (orig)', 'θ2_{hat} (orig)', 'θ3_{hat} (orig)', ...
       'θ1_{hat} (new)', 'θ2_{hat} (new)', 'θ3_{hat} (new)');
xlabel('Time (t)'); ylabel('θ̂_i'); grid on;


%% ========== SYSTEM DYNAMICS WITH LYAPUNOV SERIES PARALLEL LAW ==========
function xxdot = LyapSP(t, xx, modelType)
    global c d f gamma thetam theta_true newInput
    
    % States and estimates
    x1 = xx(1); x1_hat = xx(2);
    theta1 = xx(3); theta2 = xx(4); theta3 = xx(5);

    % Select input signal (Comment or uncomment the desired signal manually)
    if newInput
        %u = 2.5 * sin(0.05 * t);              % Testing with a lower frequency
        %u = 2.5 * sin(6 * t);                 % Testing with a higher frequency
        u = 1.5*sin(0.6*t) + 0.8*sin(6.5*t);   % Testing with a mixed frequency approach
    else
        u = c * sin(d * t) + f + 0.1 * sin(0.2*t) + 0.5 * sin(0.3*t);  % Original input
    end
    
    % Tracking error
    e1 = x1 - x1_hat;  

    % Define basis functions φ based on model type
    switch modelType
        case 1
            phi = [tanh(x1); 1/(1+x1^2); 0];
        case 2
            phi = [x1^2; sin(x1); x1^2/(x1^2+1)];
        case 3
            phi = [x1^3; tanh(x1); 1/(1+x1^2)];
        case 4
            phi = [x1; x1^2; x1^3];
        case 5
            phi = [exp(-x1^2); exp(-(x1-1)^2); exp(-(x1+1)^2)];
        otherwise
            error('Invalid model type');
    end

    % Define the system dynamics and adaptation laws
    xxdot = zeros(size(xx));
    xxdot(1) = -x1^3 + theta_true(1)*tanh(x1) + theta_true(2)/(1+x1^2) + u;     % True system
    xxdot(2) = theta1*phi(1) + theta2*phi(2) + theta3*phi(3) + u + thetam*e1;   % Estimated system
    xxdot(3) = gamma(1) * e1 * phi(1);                                          % θ1 estimate
    xxdot(4) = gamma(2) * e1 * phi(2);                                          % θ2 estimate
    xxdot(5) = gamma(3) * e1 * phi(3);                                          % θ3 estimate
end
