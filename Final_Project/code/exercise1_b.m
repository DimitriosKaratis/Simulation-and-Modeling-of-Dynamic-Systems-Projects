%% SIMULATION AND MODELING OF DYNAMIC SYSTEMS
% Final Project – Exercise 1b
% Dimitrios Karatis – 10775

clear; clc; close all;
format long

%% System Parameters
global A B c d f h i j gamma bar_omega_A bar_omega_B sigma_val M_sigma

% True system matrices
A = [-2.15  0.25;
     -0.75 -2];

B = [0;
     1.5];

% Input signal u = c * sin(d * t) + f * sin(h * t) + i * sin(j * t) coefficients
c = 5.0;
d = 0.5;
f = 4.0;
h = 1.3;
i = 3.0;
j = 2.7;

% Adaptation gains
gamma_A = 7.3;
gamma_B = 0.3;
gamma = [gamma_A, gamma_B];

% Bias disturbance upper bounds
bar_omega_A = 0.02;   
bar_omega_B = 0.05;

% Sigma modification parameters/gains
M_a11 = 2.2;
M_a12 = 0.35;
M_a21 = 0.77;
M_a22 = 2.2;
M_b1  = 0.3;
M_b2  = 2.0;

sigma_a11 = 0.03;
sigma_a12 = 0.01;
sigma_a21 = 0.02;
sigma_a22 = 0.1;
sigma_b1  = 0.02;
sigma_b2  = 0.05;


M_sigma = [M_a11 M_a12 M_a21 M_a22 M_b1 M_b2];
sigma_val = [sigma_a11 sigma_a12 sigma_a21 sigma_a22 sigma_b1 sigma_b2];

%% Simulation Settings
tStart = 0;
tStep  = 0.01;
tEnd   = 200;
tspan  = tStart:tStep:tEnd;

%% Initial Conditions
% [x1; x2; x1_hat; x2_hat; a11_hat; a12_hat; a21_hat; a22_hat; b1_hat; b2_hat; ds5; ds6; ds7; ds8; ds9; ds10]
initCond = [0; 0; 0; 0; ...
            -2; 0.2; -0.6; -1.5; ...
            -0.1; 1.0; ...
            0; 0; 0; 0; 0; 0];  

%% Run The Simulation
[t, x] = ode45(@(t, x) LyapParWithBias(t, x), tspan, initCond);

% Extract states
x1     = x(:,1);  
x2     = x(:,2);
x1_hat = x(:,3);  
x2_hat = x(:,4);

% Extract parameter estimates
a11_hat = x(:,5); 
a12_hat = x(:,6);
a21_hat = x(:,7); 
a22_hat = x(:,8);
b1_hat  = x(:,9); 
b2_hat  = x(:,10);

% Extract sigma modification parameters
sd5 = x(:,11);
sd6 = x(:,12);
sd7 = x(:,13);
sd8 = x(:,14);
sd9 = x(:,15);
sd10 = x(:,16);

% Estimation Errors
e1 = x1 - x1_hat;
e2 = x2 - x2_hat;

%% === PLOTS ===

% State Estimation
figure;
subplot(2,1,1);
plot(t, x1, 'r', t, x1_hat, 'b', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('x1');
legend('x1','x1_{hat}'); grid on;
title('State x1 vs x1_{hat}');

subplot(2,1,2);
plot(t, x2, 'r', t, x2_hat, 'b', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('x2');
legend('x2','x2_{hat}'); grid on;
title('State x2 vs x2_{hat}');

% Estimation Errors
figure;
subplot(2,1,1);
plot(t, e1, 'k', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('e1'); title('Estimation Error: e1 = x1 - x1_{hat}'); grid on;

subplot(2,1,2);
plot(t, e2, 'k', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('e2'); title('Estimation Error: e2 = x2 - x2_{hat}'); grid on;

% Parameter Estimation
figure;
subplot(3,2,1); plot(t, a11_hat, 'b'); yline(A(1,1), '--r'); title('a_{11}'); grid on;
subplot(3,2,2); plot(t, a12_hat, 'b'); yline(A(1,2), '--r'); title('a_{12}'); grid on;
subplot(3,2,3); plot(t, a21_hat, 'b'); yline(A(2,1), '--r'); title('a_{21}'); grid on;
subplot(3,2,4); plot(t, a22_hat, 'b'); yline(A(2,2), '--r'); title('a_{22}'); grid on;
subplot(3,2,5); plot(t, b1_hat, 'b');  yline(B(1), '--r');   title('b_{1}');  grid on;
subplot(3,2,6); plot(t, b2_hat, 'b');  yline(B(2), '--r');   title('b_{2}');  grid on;
sgtitle('Estimated Parameters vs True Values');

% Sigma Modification Parameter Estimation
figure;
subplot(3,2,1); plot(t, sd5, 'LineWidth', 1.2); title('σδ for a_{11}'); grid on;
subplot(3,2,2); plot(t, sd6, 'LineWidth', 1.2); title('σδ for a_{12}'); grid on;
subplot(3,2,3); plot(t, sd7, 'LineWidth', 1.2); title('σδ for a_{21}'); grid on;
subplot(3,2,4); plot(t, sd8, 'LineWidth', 1.2); title('σδ for a_{22}'); grid on;
subplot(3,2,5); plot(t, sd9, 'LineWidth', 1.2); title('σδ for b_{1}'); grid on;
subplot(3,2,6); plot(t, sd10, 'LineWidth', 1.2); title('σδ for b_{2}'); grid on;
sgtitle('Activation of σδ for Each Parameter');


%% Final Estimates Output
fprintf('\n=== Final Estimated Parameters ===\n');
fprintf('a11: True = %.4f,\t Estimated = %.4f\n', A(1,1), a11_hat(end));
fprintf('a12: True = %.4f,\t Estimated = %.4f\n', A(1,2), a12_hat(end));
fprintf('a21: True = %.4f,\t Estimated = %.4f\n', A(2,1), a21_hat(end));
fprintf('a22: True = %.4f,\t Estimated = %.4f\n', A(2,2), a22_hat(end));
fprintf('b1 : True = %.4f,\t Estimated = %.4f\n', B(1),   b1_hat(end));
fprintf('b2 : True = %.4f,\t Estimated = %.4f\n\n', B(2),   b2_hat(end));

%% Bias Variation Study
omega_levels = [0.02, 0.1, 0.3, 0.5, 1.0 5.0 10.0]; 
results = [];

for k = 1:length(omega_levels)
    bar_omega_A = omega_levels(k);  
    bar_omega_B = omega_levels(k);

    [t, x] = ode45(@(t, x) LyapParWithBias(t, x), tspan, initCond);

    e1 = x(:,1) - x(:,3);
    e2 = x(:,2) - x(:,4);

    % RMS Errors
    rms_e1 = sqrt(mean(e1.^2));
    rms_e2 = sqrt(mean(e2.^2));

    % Final parameter estimates
    final_est = x(end, 5:10);

    % Store: [omega_level, rms_e1, rms_e2, a11_hat, ..., b2_hat]
    results = [results; omega_levels(k), rms_e1, rms_e2, final_est];
end

% Display Table
disp(' ω     |   RMS(e1)   RMS(e2)   a11_hat   a12_hat   a21_hat   a22_hat   b1_hat   b2_hat');
disp('--------------------------------------------------------------------------------------');
fprintf('%5.2f |  %8.4f  %8.4f  %9.4f  %9.4f  %9.4f  %9.4f  %8.4f  %8.4f\n', results');

% === Plot RMS Estimation Errors vs Bias Level ===
figure;
plot(results(:,1), results(:,2), 'r-o', 'LineWidth', 2); hold on;
plot(results(:,1), results(:,3), 'b-o', 'LineWidth', 2);
xlabel('\omegā (Bias Upper Bound)');
ylabel('RMS Estimation Error');
legend('RMS(e_1)', 'RMS(e_2)', 'Location', 'northwest');
title('Effect of \omegā on State Estimation Error');
grid on;

% === Plot Final Parameter Estimate Errors vs Bias Level ===
param_names = {'a_{11}', 'a_{12}', 'a_{21}', 'a_{22}', 'b_{1}', 'b_{2}'};
true_params = [-2.15, 0.25, -0.75, -2.0, 0.0, 1.5];

figure;
for i = 1:6
    subplot(3,2,i);
    plot(results(:,1), results(:,i+3) - true_params(i), 'k-o', 'LineWidth', 1.5);
    yline(0, '--r');
    xlabel('\omegā'); ylabel(sprintf('Error in %s', param_names{i}));
    title(sprintf('Final Error in %s vs \\omegā', param_names{i}));
    grid on;
end
sgtitle('Parameter Estimation Errors under Increasing Bias');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xxdot = LyapParWithBias(t, xx)
    % Lyapunov Parallel Estimator with Bias and Discontinuous σ-Modification
    % Estimates all entries of A and B
    % Includes projection and robustness to bias ω(t)

    global A B c d f h i j gamma bar_omega_A bar_omega_B sigma_val M_sigma

    % System states
    x1      = xx(1);
    x2      = xx(2);
    x1_hat  = xx(3);
    x2_hat  = xx(4);

    % Parameter estimates
    a11_hat = xx(5);
    a12_hat = xx(6);
    a21_hat = xx(7);
    a22_hat = xx(8);
    b1_hat  = xx(9);
    b2_hat  = xx(10);
    
    % Input signal
    u = c * sin(d * t) + f * sin(h * t) + i * sin(j * t);

    % Bias disturbance ω(t) with bounded norm
    omega = [bar_omega_A; bar_omega_B* sin(h*t)];

    % Estimation errors
    e1 = x1 - x1_hat;
    e2 = x2 - x2_hat;

    xxdot = zeros(16,1);

    %% True system dynamics with bias
    xxdot(1) = A(1,1)*x1 + A(1,2)*x2 + B(1)*u + omega(1);
    xxdot(2) = A(2,1)*x1 + A(2,2)*x2 + B(2)*u + omega(2);

    %% Estimated system dynamics
    xxdot(3) = a11_hat * x1_hat + a12_hat * x2_hat + b1_hat * u;
    xxdot(4) = a21_hat * x1_hat + a22_hat * x2_hat + b2_hat * u;

    %% Robust Parameter adaptation with σ-modification
    % Compute σδ for each parameter
    s5  = sigma_delta(a11_hat, M_sigma(1), sigma_val(1));
    s6  = sigma_delta(a12_hat, M_sigma(2), sigma_val(2));
    s7  = sigma_delta(a21_hat, M_sigma(3), sigma_val(3));
    s8  = sigma_delta(a22_hat, M_sigma(4), sigma_val(4));
    s9  = sigma_delta(b1_hat,  M_sigma(5), sigma_val(5));
    s10 = sigma_delta(b2_hat, M_sigma(6), sigma_val(6));

    % Gradient + σ-modification adaptation
    xxdot(5)  = gamma(1) * e1 * x1_hat - gamma(1) * s5  * a11_hat;
    xxdot(6)  = gamma(1) * e1 * x2_hat - gamma(1) * s6  * a12_hat;
    xxdot(7)  = gamma(1) * e2 * x1_hat - gamma(1) * s7  * a21_hat;
    xxdot(8)  = gamma(1) * e2 * x2_hat - gamma(1) * s8  * a22_hat;
    xxdot(9)  = gamma(2) * e1 * u      - gamma(2) * s9  * b1_hat;
    xxdot(10) = gamma(2) * e2 * u      - gamma(2) * s10 * b2_hat;

    % Store delta_sigma values into the state vector for plotting
    xxdot(11) = s5;
    xxdot(12) = s6;
    xxdot(13) = s7;
    xxdot(14) = s8;
    xxdot(15) = s9;
    xxdot(16) = s10;


    %% Constraint projection
    % Constraint: a11 ∈ [-3, -1]
    if a11_hat < -3 || a11_hat > -1
        xxdot(5) = 0;
    elseif a11_hat == -3 && xxdot(5) < 0
        xxdot(5) = 0;
    elseif a11_hat == -1 && xxdot(5) > 0
        xxdot(5) = 0;
    end

    % Constraint: b2 ≥ 1
    if b2_hat < 1
        xxdot(10) = 0;
    elseif b2_hat == 1 && xxdot(10) < 0
        xxdot(10) = 0;
    end
end

%% Helper function: Continuous σ-modification term
function sigma_d = sigma_delta(theta_hat, M, sigma)
    abs_theta = abs(theta_hat);
    if abs_theta < M
        sigma_d = 0;
    elseif (abs_theta <= 2*M && abs_theta >= M)
        sigma_d = sigma * (abs_theta / M - 1);
    elseif (abs_theta > 2*M)
        sigma_d = sigma;
    end
end
