%% SIMULATION AND MODELING OF DYNAMIC SYSTEMS
% Assignment 2 - April 2025
% KARATIS DIMITRIOS 10775

%% Exercise 1 - a
clear; clc; close all;

% Global variables
global gamma m b k c d f p1 p2

% True system parameters
m = 1.315;         % Mass
b = 0.225;         % Damping
k = 0.725;         % Stiffness

% Filter and adaptation gains
gamma = 0.05;
lambda1 = 1.0;
lambda2 = lambda1;

% Solve the polynomial equation p1^2 - lambda1*p1 + lambda1 = 0
coeffs = [1 -lambda1 lambda1];  
sol = roots(coeffs);

% Check sol(1)
if real(sol(1)) > 0 && real(lambda1 - sol(1)) > 0
    p1 = sol(1);
    p2 = lambda1 - p1;
    
% If not, check sol(2)
elseif real(sol(2)) > 0 && real(lambda1 - sol(2)) > 0
    p1 = sol(2);
    p2 = lambda1 - p1;

% If neither work, throw an error
else
    error('No valid root found where both p1 and p2 have positive real parts.');
end

% Simulation time
tspan = 0:0.01:20;

% Initial conditions: [x; x_dot; j1; j2; j3; dj1; dj2; dj3; theta1; theta2; theta3]
initCond = [0; 0; 0; 0; 0; 0; 0; 0; 0.01; 0.01; 0.01];

% Loop over different sets of input parameters
for input_set = 1:2
    if input_set == 1
        % Set parameters for u = 2.5 -> u = csin(dt) + f
        c = 0;         
        d = 0;           
        f = 2.5;       
    else
        % Set parameters for u = 2.5sin(t) -> u = csin(dt) + f
        c = 2.5;          
        d = 1;            
        f = 0;         
    end

    % Solve the ODE
    [t, x] = ode45(@(t, x) GraDE_2nd_order(t, x), tspan, initCond);

    % Extract variables
    x_sketo = x(:,1);
    x_dot= x(:, 2);
    j_est1 = x(:, 3);
    j_est2 = x(:, 4);
    j_est3 = x(:, 5);
    dj_est1 = x(:, 6);
    dj_est2 = x(:, 7);
    dj_est3 = x(:, 8);
    theta_est1 = x(:, 9);
    theta_est2 = x(:, 10);
    theta_est3 = x(:, 11);

    theta_est = [theta_est1, theta_est2, theta_est3];
    j_est = [j_est1, j_est2, j_est3];

    % True parameter vector
    theta_truee = [b/m - lambda1, k/m - lambda2, 1/m];

    % Estimation of x
    x_est = theta_est1.*j_est1 + theta_est2.*j_est2 + theta_est3.*j_est3;

    % X error
    e_x = x_sketo - x_est;

    m_est = 1 / theta_est(end, 3);
    k_est = (theta_est(end,2) + lambda2) * m_est;
    b_est = (theta_est(end, 1) + lambda1) * m_est;

    % Plot x(t) and x_hat(t)
    figure;
    plot(t, x_sketo, 'b', 'LineWidth', 1.5); hold on;
    plot(t, x_est, 'r--', 'LineWidth', 1.5);
    legend('x(t)', 'x̂(t)'); xlabel('Time [s]'); ylabel('Output'); grid on;
    title(['Actual vs Estimated Output for input set ', num2str(input_set)]);

    % Plot x error
    figure;
    plot(t, e_x, 'k', 'LineWidth', 1.5);
    xlabel('Time [s]'); ylabel('e_x(t) = x(t) - x̂(t)'); grid on;
    title(['Estimation Error for input set ', num2str(input_set)]);

    m_hat = 1 ./ theta_est(:,3);
    b_hat = (theta_est(:,1) + lambda1) .* m_hat;
    k_hat = (theta_est(:,2) + lambda2) .* m_hat;

    % Mask for t >= 0.5
    mask = t >= 0.5;

    % Apply mask to all relevant vectors
    t_plot = t(mask);
    m_hat_plot = m_hat(mask);
    b_hat_plot = b_hat(mask);
    k_hat_plot = k_hat(mask);

    % Plot m_est and true m
    figure;
    subplot(3,1,1);
    plot(t, m_hat, 'LineWidth', 1.5); 
    yline(m, '--r', 'LineWidth', 1.5); 
    ylabel('m_est(t)'); 
    grid on; 
    ylim([1.5 2.5]);

    % Plot b_est and true b
    subplot(3,1,2);
    plot(t, b_hat, 'LineWidth', 1.5); 
    yline(b, '--r', 'LineWidth', 1.5); 
    ylabel('b_est(t)'); 
    grid on;

    % Plot k_est and true k
    subplot(3,1,3);
    plot(t, k_hat, 'LineWidth', 1.5); 
    yline(k, '--r', 'LineWidth', 1.5);  
    ylabel('k_est(t)'); 
    xlabel('Time [s]'); 
    grid on;
    sgtitle(['Estimated Parameters over Time for input set ', num2str(input_set)]);

     % ==== Final estimated values ====
    if input_set == 1
        fprintf('\n=== Final Parameter Estimates for Set 1 (u = 2.5) ===\n');
    else
        fprintf('\n=== Final Parameter Estimates for Set 2 (u = 2.5sin(t)) ===\n');
    end
    fprintf('Mass:      m_true = %.4f\t m_est = %.4f\n', m, m_est);
    fprintf('Damping:   b_true = %.4f\t b_est = %.4f\n', b, b_est);
    fprintf('Stiffness: k_true = %.4f\t k_est = %.4f\n', k, k_est);
end




