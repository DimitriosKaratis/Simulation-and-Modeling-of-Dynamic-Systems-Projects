%% SIMULATION AND MODELING OF DYNAMIC SYSTEMS
% Assignment 1 - March 2025
% KARATIS DIMITRIOS 10775

%% Exercise 2 - a

clear all;
close all;
clc;

% Define the parameters of the system
m = 0.75;  % kg
L = 1.25; % m
c = 0.15; % N*m*sec
g = 9.81; % m/sec^2
A0 = 4;   % Ν*m
omega = 2; % rad/sec

%% Simulation of the REAL system
% Define the simulation time
T_sim = 20; 
dt = 1e-4;  
tspan = 0:dt:T_sim;

% Initial conditions
x0 = [0; 0];

% Input function
u = @(t) A0 * sin(omega * t);

% Define the system
dynamics = @(t, x) [x(2); (1/(m*L^2)) * (u(t) - c*x(2) - m*g*L*x(1))];

% Solve the system using ode45
[t, X] = ode45(dynamics, tspan, x0);

q = X(:,1);
q_dot = X(:,2);

%% Estimating with Least-Squares method

% Pole of the system
p = 0.5;

% Downsampling to 0.1-second intervals
dt_sampled = 0.1;  
t_sampled = 0:dt_sampled:T_sim;  

% Interpolate the high-resolution data to match the new sampling points
q_sampled = interp1(t, q, t_sampled);  
q_dot_sampled = interp1(t, q_dot, t_sampled);  
u_sampled = A0 * sin(omega * t_sampled); 

% Creating the phi matrix
phi = zeros(length(t_sampled'),3);

phi1 = lsim(tf(-1,[1 p]), q_dot_sampled', t_sampled');
phi2 = lsim(tf(-1,[1 p]), q_sampled', t_sampled');
phi3 = lsim(tf(1,[1 p]), u_sampled', t_sampled');   

phi(:,1) = phi1;
phi(:,2) = phi2;
phi(:,3) = phi3;

% Phi squared (or ΦΤΦ)
phiTphi = (phi') * phi;                              

theta0 = inv(phiTphi) * (phi') * (q_dot_sampled');

% Find the estimated parameters
L_estim = g / (theta0(2));
m_estim = 1 / ((L_estim^2) * (theta0(3)));
c_estim = (m_estim * L_estim^2) * (theta0(1) + p);

% Display the original and estimated values
fprintf('\n--- Original System Parameters ---\n');
fprintf('Mass (m): %.2f kg\n', m);
fprintf('Length (L): %.2f m\n', L);
fprintf('Damping coefficient (c): %.2f N*m*sec\n', c);

fprintf('\n--- Estimated System Parameters (Least-Squares) ---\n');
fprintf('Estimated Mass (m_estim): %.2f kg\n', m_estim);
fprintf('Estimated Length (L_estim): %.2f m\n', L_estim);
fprintf('Estimated Damping coefficient (c_estim): %.2f N*m*sec\n', c_estim);

%% Simulate the new response of the system based on the newly estimated parameters

% Define the parameters of the system
m = m_estim;  % kg
L = L_estim; % m
c = c_estim; % N*m*sec
g = 9.81; % m/sec^2
A0 = 4;   % Ν*m
omega = 2; % rad/sec

% Define the simulation time
T_sim = 20; 
dt = 1e-4; 
tspan = 0:dt:T_sim;

% Initial conditions
x0 = [0; 0];

% Input function
u = @(t) A0 * sin(omega * t);

% Define the system
dynamics = @(t, x) [x(2); (1/(m*L^2)) * (u(t) - c*x(2) - m*g*L*x(1))];

% Solve the system using ode45
[t, X] = ode45(dynamics, tspan, x0);

q_estim = X(:,1);
q_dot_estim = X(:,2);

% Plots
figure;
subplot(2,1,1);
plot(t, q_estim, 'b', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('q [rad]'); grid on;
title('Angular displacement response (WITH ESTIMATED PARAMETERS)');

subplot(2,1,2);
plot(t, q_dot_estim, 'r', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('dq/dt [rad/s]'); grid on;
title('Angular velocity response (WITH ESTIMATED PARAMETERS)');

%% Plot the diagramms
% Calculate error between original and estimated responses
error_q = q - q_estim;
mse_error = mean(error_q.^2);

fprintf('\n\nMean Squared Error : %.8f\n', mse_error);

% Create the combined plot
figure;
hold on;

% Original response
plot(t, q, 'b-', 'LineWidth', 1.5);

% Estimated response
plot(t, q_estim, 'r--', 'LineWidth', 1.5);

% Error plot
plot(t, error_q, 'g-.', 'LineWidth', 1.5);

% Add legends and labels
legend('Original q(t)', 'Estimated q_{estimate}(t)', 'Error e(t) = q(t) - q_{estimate}(t)');
xlabel('Time [s]');
ylabel('Response [rad]');
title('Original vs Estimated Response and Error');
grid on;

hold off;

