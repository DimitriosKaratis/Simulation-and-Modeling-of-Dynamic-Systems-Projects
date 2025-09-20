%% SIMULATION AND MODELING OF DYNAMIC SYSTEMS
% Assignment 2 - April 2025
% KARATIS DIMITRIOS 10775

%% Exercise 2 - b
%% Lyapunov Method / Real time estimator / Series-Parallel Configuration 

% Clear environment
clear; clc; close all;

% Time settings
dt = 1e-4;          % Integration step
time = 0:dt:20;     % Time vector

% Desired reference signal (command input)
T_final = 20;
r_d = @(t) (pi/20)*(1 - cos(2*pi*t/T_final));  % Desired trajectory

% System actual parameters (ground truth)
true_a1 = 1.315;
true_a2 = 0.725;
true_a3 = 0.225;
true_b  = 1.175;

% Control design parameters
phi_initial = 1.0;
phi_final = 0.01;
decay_rate = 1;
rho_val = 1.5;
k1_val = 2;
k2_val = 5;

% Initial states [r, r_dot]
x0 = [0, 0];

% Simulate system response
[t_vec, x_sim] = ode45(@(t, x) roll_angle_dynamics(t, x, true_a1, true_a2, true_a3, true_b, ...
    phi_initial, phi_final, decay_rate, rho_val, k1_val, k2_val, r_d), time, x0);

% Extract simulation data
r_vals = x_sim(:, 1);
r_dot_vals = x_sim(:, 2);
r_d_vals = arrayfun(r_d, t_vec);

% Compute control input over time
u_input = zeros(1, length(r_vals));
for idx = 1:length(r_vals)
    t_now = t_vec(idx);
    r_now = r_vals(idx);
    r_dot_now = r_dot_vals(idx);
    r_desired = r_d_vals(idx);

    phi_t = (phi_initial - phi_final) * exp(-decay_rate * t_now) + phi_final;
    z1 = (r_now - r_desired) / phi_t;
    a_val = -k1_val * log((1 + z1) / (1 - z1));

    z2 = (r_dot_now - a_val) / rho_val;
    u_input(idx) = -k2_val * log((1 + z2) / (1 - z2));
end

%% Parameter Estimation via Series-Parallel Approach

% Learning rates
gamma1 = 8.5;
gamma2 = 0.027;
gamma3 = 1.15;
gamma4 = 0.03;

% Observer gains
theta1 = 50;
theta2 = 90;

% Initialize parameter estimates
est_a1 = zeros(1, length(r_vals));
est_a2 = zeros(1, length(r_vals));
est_a3 = zeros(1, length(r_vals));
est_b  = zeros(1, length(r_vals));

est_r = zeros(1, length(r_vals));
est_r(1) = r_vals(1);

est_r_dot = zeros(1, length(r_vals));
est_r_dot(1) = r_dot_vals(1);

% Initial guesses
est_a1(1) = 1.0;
est_a2(1) = 1.0;
est_a3(1) = 0.5;
est_b(1)  = 1.0;

% Adaptation loop
for k = 1:length(r_vals)-1
    err_dot = r_dot_vals(k) - est_r_dot(k);

    f_func = sin(r_vals(k));
    g_func = r_dot_vals(k)^2 * sin(2 * r_vals(k));

    % Update laws (steepest descent)
    da1 = -gamma1 * r_dot_vals(k) * err_dot;
    da2 = -gamma2 * f_func * err_dot;
    da3 =  gamma3 * g_func * err_dot;
    db  =  gamma4 * u_input(k) * err_dot;

    est_a1(k+1) = est_a1(k) + da1;
    est_a2(k+1) = est_a2(k) + da2;
    est_a3(k+1) = est_a3(k) + da3;
    est_b(k+1)  = est_b(k)  + db;

    % Observer dynamics
    est_r_ddot = - est_a1(k) * r_dot_vals(k) - est_a2(k) * f_func + est_a3(k) * g_func + est_b(k) * u_input(k) ...
                 + theta2 * (r_dot_vals(k) - est_r_dot(k));

    est_r_dot(k+1) = est_r_dot(k) + dt * est_r_ddot;
    est_r(k+1) = est_r(k) + dt * (est_r_dot(k) + theta1 * (r_vals(k) - est_r(k)));
end

% Estimation error
r_error = r_vals - est_r';

%% Plotting Results

% True vs Estimated displacement
figure;
plot(t_vec, r_vals, 'b', 'LineWidth', 1.5);
hold on;
plot(t_vec, est_r, 'r--', 'LineWidth', 1.5);
legend('True r', 'Estimated r_est', 'Location', 'best');
xlabel('Time [s]');
ylabel('Roll Angle [rad]');
title('Roll Angle: True vs Estimated');

% Estimation Error
figure;
plot(t_vec, r_error, 'k', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Estimation Error [rad]');
title('Estimation Error: r(t) - r_est(t)');

% Parameter estimation over time
figure;
subplot(4,1,1);
plot(t_vec, est_a1, 'LineWidth', 1.5); hold on;
yline(true_a1, '--k', 'True a1');
legend('Estimated a1', 'Location', 'best');
ylabel('a1');
xlabel('Time [s]');
title('Estimated a1');

subplot(4,1,2);
plot(t_vec, est_a2, 'LineWidth', 1.5); hold on;
yline(true_a2, '--k', 'True a2');
legend('Estimated a2', 'Location', 'best');
ylabel('a2');
xlabel('Time [s]');
title('Estimated a2');

subplot(4,1,3);
plot(t_vec, est_a3, 'LineWidth', 1.5); hold on;
yline(true_a3, '--k', 'True a3');
legend('Estimated a3', 'Location', 'best');
ylabel('a3');
xlabel('Time [s]');
title('Estimated a3');

subplot(4,1,4);
plot(t_vec, est_b, 'LineWidth', 1.5); hold on;
yline(true_b, '--k', 'True b');
legend('Estimated b', 'Location', 'best');
ylabel('b');
xlabel('Time [s]');
title('Estimated b');

%% Final Estimated Values
fprintf('\n=== Final Parameter Estimates (Series-Parallel Configuration) ===\n');
fprintf('a1_true = %.4f\t a1_est = %.4f\n', true_a1, est_a1(end));
fprintf('a2_true = %.4f\t a2_est = %.4f\n', true_a2, est_a2(end));
fprintf('a3_true = %.4f\t a3_est = %.4f\n', true_a3, est_a3(end));
fprintf('b_true = %.4f\t         b_est = %.4f\n', true_b, est_b(end));