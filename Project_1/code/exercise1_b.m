%% SIMULATION AND MODELING OF DYNAMIC SYSTEMS
% Assignment 1 - March 2025
% KARATIS DIMITRIOS 10775

%% Exercise 1 - b

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

% Plot the response of the system
transfer_f = tf(1, [m*(L^2) c g*m*L]);
lsim(transfer_f, u(t), t);

% Plots
figure;
subplot(2,1,1);
plot(t, q, 'b', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('q [rad]'); grid on;
title('Angular displacement response');
 
subplot(2,1,2);
plot(t, q_dot, 'r', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('dq/dt [rad/s]'); grid on;
title('Angular velocity response');


