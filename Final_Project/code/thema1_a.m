%% SIMULATION AND MODELING OF DYNAMIC SYSTEMS
% Final Project – Exercise 1a
% Dimitrios Karatis – 10775

clear; clc; close all;
format long

%% System Parameters
global A B c d f h i j gamma

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
gamma_A = 16;
gamma_B = 3;
gamma = [gamma_A, gamma_B];

%% Simulation parameters
tStart = 0;
tStep  = 0.01;
tEnd   = 500;
tspan  = tStart:tStep:tEnd;

%% Initial Conditions
% [x1; x2; x1_hat; x2_hat; a11_hat; a12_hat; a21_hat; a22_hat; b1_hat; b2_hat]
initCond = [0; 0; 0; 0; ...
            -2; 0.2; -0.6; -1.5; ...
            -0.1; 1.0];  

%% Run simulation
[t, x] = ode45(@(t, x) LyapPar(t, x), tspan, initCond);

% Extract system states
x1 = x(:,1);     
x2 = x(:,2);
x1_hat = x(:,3); 
x2_hat = x(:,4);

% Extract parameter estimates
a11_hat = x(:,5); 
a12_hat = x(:,6);
a21_hat = x(:,7); 
a22_hat = x(:,8);
b1_hat  = x(:,9); 
b2_hat  = x(:,10);

% Compute estimation errors
e1 = x1 - x1_hat;
e2 = x2 - x2_hat;

%% === PLOTS ===

%% State Estimation Plot
figure;
subplot(2,1,1);
plot(t, x1, 'r', t, x1_hat, 'b', 'LineWidth', 1.5);
title('x1 vs x1_{hat} (Parallel Estimator)');
xlabel('Time [s]'); ylabel('State'); grid on;
legend('x1','x1_{hat}');

subplot(2,1,2);
plot(t, x2, 'r', t, x2_hat, 'b', 'LineWidth', 1.5);
title('x2 vs x2_{hat} (Parallel Estimator)');
xlabel('Time [s]'); ylabel('State'); grid on;
legend('x_2','x2_{hat}');

%% Estimation Error Plot
figure;
subplot(2,1,1);
plot(t, e1, 'k', 'LineWidth', 1.5);
title('Estimation Error e1 = x1 - x1_{hat}');
xlabel('Time [s]'); ylabel('Error'); grid on;

subplot(2,1,2);
plot(t, e2, 'k', 'LineWidth', 1.5);
title('Estimation Error e2 = x2 - x2_{hat}');
xlabel('Time [s]'); ylabel('Error'); grid on;

%% Parameter Estimation Plot
figure;
subplot(3,2,1); plot(t, a11_hat, 'b'); yline(A(1,1),'--r'); title('a_{11}'); grid on;
subplot(3,2,2); plot(t, a12_hat, 'b'); yline(A(1,2),'--r'); title('a_{12}'); grid on;
subplot(3,2,3); plot(t, a21_hat, 'b'); yline(A(2,1),'--r'); title('a_{21}'); grid on;
subplot(3,2,4); plot(t, a22_hat, 'b'); yline(A(2,2),'--r'); title('a_{22}'); grid on;
subplot(3,2,5); plot(t, b1_hat, 'b');  yline(B(1),'--r');   title('b_{1}');  grid on;
subplot(3,2,6); plot(t, b2_hat, 'b');  yline(B(2),'--r');   title('b_{2}');  grid on;

sgtitle('Parameter Estimations vs True Values');

%% Final Estimates Output
fprintf('\n=== Final Estimated Parameters ===\n');
fprintf('a11: True = %.4f,\t Estimated = %.4f\n', A(1,1), a11_hat(end));
fprintf('a12: True = %.4f,\t Estimated = %.4f\n', A(1,2), a12_hat(end));
fprintf('a21: True = %.4f,\t Estimated = %.4f\n', A(2,1), a21_hat(end));
fprintf('a22: True = %.4f,\t Estimated = %.4f\n', A(2,2), a22_hat(end));
fprintf('b1 : True = %.4f,\t Estimated = %.4f\n', B(1),   b1_hat(end));
fprintf('b2 : True = %.4f,\t Estimated = %.4f\n', B(2),   b2_hat(end));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xxdot = LyapPar(t, xx)
    % Lyapunov Parallel Estimator for Exercise 1a
    % Estimates all entries of A and B
    % with projection constraints on a11 and b2

    global A B c d f h i j gamma

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

    % Estimation errors
    e1 = x1 - x1_hat;
    e2 = x2 - x2_hat;

    xxdot = zeros(10,1);

    %% True system dynamics
    xxdot(1) = A(1,1)*x1 + A(1,2)*x2 + B(1)*u;
    xxdot(2) = A(2,1)*x1 + A(2,2)*x2 + B(2)*u;

    %% Estimated system dynamics
    xxdot(3) = a11_hat * x1_hat + a12_hat * x2_hat + b1_hat * u;
    xxdot(4) = a21_hat * x1_hat + a22_hat * x2_hat + b2_hat * u;

    %% Parameter adaptation laws (gradient update)
    xxdot(5)  = gamma(1) * e1 * x1_hat;   % a11_hat_dot
    xxdot(6)  = gamma(1) * e1 * x2_hat;   % a12_hat_dot
    xxdot(7)  = gamma(1) * e2 * x1_hat;   % a21_hat_dot
    xxdot(8)  = gamma(1) * e2 * x2_hat;   % a22_hat_dot
    xxdot(9)  = gamma(2) * e1 * u;        % b1_hat_dot
    xxdot(10) = gamma(2) * e2 * u;        % b2_hat_dot

    %% Constraint projection

    % Constraint: a11 in [-3, -1]
    if a11_hat < -3 || a11_hat > -1
        xxdot(5) = 0;
    elseif a11_hat == -3 && xxdot(5) < 0
        xxdot(5) = 0;
    elseif a11_hat == -1 && xxdot(5) > 0
        xxdot(5) = 0;
    end

    % Constraint: b2 >= 1
    if b2_hat < 1
        xxdot(10) = 0;
    elseif b2_hat == 1 && xxdot(10) < 0
        xxdot(10) = 0;
    end

  

end
