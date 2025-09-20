%% SIMULATION AND MODELING OF DYNAMIC SYSTEMS
% Assignment 1 - March 2025
% KARATIS DIMITRIOS 10775

%% Exercise 3 - a

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

% Simulation setup
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

% Downsampling setup
dt_sampled = 0.1;  
t_sampled = 0:dt_sampled:T_sim;  
q_sampled = interp1(t, q, t_sampled);
u_sampled = A0 * sin(omega * t_sampled);

% Custom noise parameters
custom_mean = 0;
custom_std_dev = 0.1;

% Filter setup
p1 = 0.5;
p2 = 0.5;

% Preallocate vectors for estimated parameters
num_trials = 100;
L_estimations = zeros(1, num_trials);
m_estimations = zeros(1, num_trials);
c_estimations = zeros(1, num_trials);

for i = 1:num_trials
    % Add Gaussian noise
    q_noisy = q_sampled + normrnd(custom_mean, custom_std_dev, size(q_sampled));
    u_noisy = u_sampled + normrnd(custom_mean, custom_std_dev, size(u_sampled));

    % Construct phi matrix
    phi_noisy = zeros(length(t_sampled'), 3);
    phi_noisy(:,1) = lsim(tf([-1 0],[1 (p1+p2) p1*p2]), q_noisy', t_sampled');
    phi_noisy(:,2) = lsim(tf(-1,[1 (p1+p2) p1*p2]), q_noisy', t_sampled');
    phi_noisy(:,3) = lsim(tf(1,[1 (p1+p2) p1*p2]), u_noisy', t_sampled');

    % Calculate theta
    phiTphi_noisy = (phi_noisy') * phi_noisy;                              
    theta0_noisy = inv(phiTphi_noisy) * (phi_noisy') * (q_noisy');

    % Estimate parameters
    L_estimations(i) = g / (theta0_noisy(2) + (p1 * p2));
    m_estimations(i) = 1 / ((L_estimations(i)^2) * (theta0_noisy(3)));
    c_estimations(i) = (m_estimations(i) * L_estimations(i)^2) * (theta0_noisy(1) + (p1 + p2));
end

% Compute mean values of the estimated parameters
L_mean = mean(L_estimations);
m_mean = mean(m_estimations);
c_mean = mean(c_estimations);

% Display results
fprintf('Estimated L (mean): %.4f m\n', L_mean);
fprintf('Real L: %.4f m\n', L);
fprintf('Estimated m (mean): %.4f kg\n', m_mean);
fprintf('Real m: %.4f kg\n', m);
fprintf('Estimated c (mean): %.4f N*m*sec\n', c_mean);
fprintf('Real c: %.4f N*m*sec\n', c);