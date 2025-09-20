%% SIMULATION AND MODELING OF DYNAMIC SYSTEMS
% Assignment 2 - April 2025
% KARATIS DIMITRIOS 10775

%% Exercise 2 - c
%% Lyapunov Method / Real-time Estimator / Series-Parallel Configuration with NOISE

% Clear environment
clear; clc; close all;

%% Time and Input Signal Setup
dt = 1e-4;                  % Integration time step
time = 0:dt:20;             % Time vector
T_final = 20;                 
r_d = @(t) (pi/20)*(1 - cos(2*pi*t/T_final));  % Desired trajectory

% Disturbance signal
disturbance = 0.15 * sin(0.5 * time);

%% True System Parameters
a1_true = 1.315;
a2_true = 0.725;
a3_true = 0.225;
b_true  = 1.175;

%% Controller Design Parameters
phi_0 = 1.0;
phi_inf = 0.01;
lambda = 1;
rho = 1.5;
k1 = 2;
k2 = 5;

% Initial conditions [r, r_dot]
x0 = [0, 0];

%% Simulate Nonlinear System Dynamics with Disturbance
[t_vec, x] = ode45(@(t, x) roll_angle_dynamics_noise(t, x, a1_true, a2_true, a3_true, b_true, ...
    phi_0, phi_inf, lambda, rho, k1, k2, r_d, disturbance, time), time, x0);

% Extract system states
r_vals = x(:, 1);
r_dot_vals = x(:, 2);
r_d_vals = arrayfun(r_d, t_vec);

%% Compute Control Input Over Time
u_vals = zeros(1, length(t_vec));
for i = 1:length(t_vec)
    t_i = t_vec(i);
    x1 = r_vals(i);
    x2 = r_dot_vals(i);
    r_des = r_d_vals(i);

    phi_t = (phi_0 - phi_inf) * exp(-lambda * t_i) + phi_inf;
    z1 = (x1 - r_des) / phi_t;
    a_val = -k1 * log((1 + z1) / (1 - z1));

    z2 = (x2 - a_val) / rho;
    u_vals(i) = -k2 * log((1 + z2) / (1 - z2));
end

%% Real-time Parameter Estimation - Series Parallel Configuration
% Learning rates
gamma1 = 8.5;
gamma2 = 0.027;
gamma3 = 1.15;
gamma4 = 0.03;

% Observer gains
theta1 = 50;
theta2 = 90;

% Initialization
a1_est = zeros(1, length(t_vec));
a2_est = zeros(1, length(t_vec));
a3_est = zeros(1, length(t_vec));
b_est  = zeros(1, length(t_vec));
r_hat  = zeros(1, length(t_vec));
r_dot_hat = zeros(1, length(t_vec));
r_hat(1) = r_vals(1);
r_dot_hat(1) = r_dot_vals(1);

% Initial parameter estimates
a1_est(1) = 1.0;
a2_est(1) = 1.0;
a3_est(1) = 0.5;
b_est(1)  = 1.0;

% Estimation loop
for i = 1:length(t_vec)-1
    e_dot = r_dot_vals(i) - r_dot_hat(i);
    f = sin(r_vals(i));
    g = r_dot_vals(i)^2 * sin(2 * r_vals(i));

    a1_dot = -gamma1 * r_dot_vals(i) * e_dot;
    a2_dot = -gamma2 * f * e_dot;
    a3_dot = gamma3 * g * e_dot;
    b_dot  = gamma4 * u_vals(i) * e_dot;

    a1_est(i+1) = a1_est(i) + a1_dot;
    a2_est(i+1) = a2_est(i) + a2_dot;
    a3_est(i+1) = a3_est(i) + a3_dot;
    b_est(i+1)  = b_est(i)  + b_dot;

    r_ddot_hat = -a1_est(i) * r_dot_vals(i) - a2_est(i) * f + a3_est(i) * g + ...
                  b_est(i) * u_vals(i) + theta2 * (r_dot_vals(i) - r_dot_hat(i)) + disturbance(i);

    r_dot_hat(i+1) = r_dot_hat(i) + dt * r_ddot_hat;
    r_hat(i+1) = r_hat(i) + dt * (r_dot_hat(i) + theta1 * (r_vals(i) - r_hat(i)));
end

% Estimation error
r_error = r_vals - r_hat';

%% Plot Results

% True vs Estimated Displacement
figure;
plot(t_vec, r_vals, 'b', 'LineWidth', 1.5); hold on;
plot(t_vec, r_hat, 'r--', 'LineWidth', 1.5);
legend('True r', 'Estimated r_est', 'Location', 'best');
xlabel('Time [s]');
ylabel('Roll Angle [rad]');
title('Roll Angle: True vs Estimated (NOISE)');

% Estimation Error
figure;
plot(t_vec, r_error, 'k', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Estimation Error [rad]');
title('Estimation Error: r(t) - r_est(t) (NOISE)');

% Parameter Estimates
figure;
subplot(4,1,1)
plot(t_vec, a1_est, 'LineWidth', 1.5); hold on; yline(a1_true, '--k');
yline(a1_true, '--k', 'True a1');
legend('Estimated a1', 'Location', 'best');
ylabel('a1');
xlabel('Time [s]');
title('Estimated a1 (NOISE)');

subplot(4,1,2)
plot(t_vec, a2_est, 'LineWidth', 1.5); hold on; yline(a2_true, '--k');
yline(a2_true, '--k', 'True a2');
legend('Estimated a2', 'Location', 'best');
ylabel('a2');
xlabel('Time [s]');
title('Estimated a2 (NOISE)');

subplot(4,1,3)
plot(t_vec, a3_est, 'LineWidth', 1.5); hold on; yline(a3_true, '--k');
yline(a3_true, '--k', 'True a3');
legend('Estimated a3 (NOISE)', 'Location', 'best');
ylabel('a3');
xlabel('Time [s]');
title('Estimated a3 (NOISE)');

subplot(4,1,4)
plot(t_vec, b_est, 'LineWidth', 1.5); hold on; yline(b_true, '--k');
yline(b_true, '--k', 'True b');
legend('Estimated b', 'Location', 'best');
ylabel('b');
xlabel('Time [s]');
title('Estimated b (NOISE)');


%% Final Estimated Values
fprintf('\n=== Final Parameter Estimates (Series-Parallel Configuration) (NOISE) ===\n');
fprintf('a1_true = %.4f\t a1_est = %.4f\n', a1_true, a1_est(end));
fprintf('a2_true = %.4f\t a2_est = %.4f\n', a2_true, a2_est(end));
fprintf('a3_true = %.4f\t a3_est = %.4f\n', a3_true, a3_est(end));
fprintf('b_true = %.4f\t         b_est = %.4f\n', b_true, b_est(end));


