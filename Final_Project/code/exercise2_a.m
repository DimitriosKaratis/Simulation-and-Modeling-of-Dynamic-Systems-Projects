%% SIMULATION AND MODELING OF DYNAMIC SYSTEMS
% Final Project – Exercise 2a
% Dimitrios Karatis – 10775

clear; clc; close all;
format long

global c d f gamma thetam theta_true

% True system dynamics
% The true function is: f(x, u, θ) = -x^3 + θ1 * tanh(x) + θ2 / (1 + x^2) + u
theta_true = [0.8, 1.2]; 

% Estimator gains
thetam = 41.2;                  
gamma = [39.9 41.7 43.3];       

% Simulation settings
tStart = 0; tEnd = 200; tStep = 0.01;
tspan = tStart:tStep:tEnd;     
N = length(tspan);            

% Input signal u = c * sin(d * t) + f + 0.1 * sin(0.2*t) + 0.5 * sin(0.3*t) parameters
c = 2.5; d = 1; f = 0;

% Energy of tracking error for each model
error_energy = zeros(1,5);      

for modelType = 1:5
    fprintf('\n===== Running for Model %d =====\n', modelType);
    
    % Initial conditions: [x1, x1_hat, θ1_hat, θ2_hat, θ3_hat]
    initCond = [0; 0; 0.01; 0.01; 0.01];
    
    % Solve the ODE system using the Lyapunov-Series-Parallel adaptation law
    [t, x] = ode45(@(tpar, xpar) LyapSP(tpar, xpar, modelType), tspan, initCond);
    
    % Extract system states and parameter estimates
    x1 = x(:,1); x1_hat = x(:,2);
    theta1_hat = x(:,3); theta2_hat = x(:,4); theta3_hat = x(:,5);
    
    % Compute tracking error
    e1 = x1 - x1_hat;
    
    % Compute energy of the tracking error (used to evaluate performance)
    error_energy(modelType) = trapz(t, e1.^2);
    
    
    figure('Name', sprintf('Model %d Results', modelType), 'NumberTitle', 'off');
    
    % Subplot 1: x1 and x1_hat
    subplot(3,1,1)
    plot(t, x1, 'r', t, x1_hat, 'b'); grid on;
    title(sprintf('State Tracking for Model %d', modelType));
    legend('x_1', 'x1_{hat}');
    xlabel('Time (t)'); ylabel('x_1');
    
    % Subplot 2: Tracking Error
    subplot(3,1,2)
    plot(t, e1, 'k'); grid on;
    title('Tracking Error e_1(t)');
    xlabel('Time (t)'); ylabel('Error');
    
    % Subplot 3: Parameter Estimates
    subplot(3,1,3)
    plot(t, theta1_hat, 'LineWidth', 1.2); hold on;
    plot(t, theta2_hat, 'LineWidth', 1.2);
    plot(t, theta3_hat, 'LineWidth', 1.2); grid on;
    title('Parameter Estimates θi_{hat}(t)');
    xlabel('Time (t)'); ylabel('θi_{hat}');
    legend('θ1_{hat}','θ2_{hat}','θ3_{hat}');

end

fprintf('\n===== Model Comparison (Full Data) =====\n');
for i = 1:5
    fprintf('Model %d: Energy = %.5f\n', i, error_energy(i));
end
[~, best] = min(error_energy);
fprintf('\nBest model: Model %d\n', best);

%% ========== M-FOLD CROSS VALIDATION ==========
M = 10;                           % Number of folds
cv_energy = zeros(1,5);           % Average validation error per model
fprintf('\n\n===== Cross Validation (%d-fold) =====\n', M);

for modelType = 1:5
    fold_energy = zeros(1,M);     

    for k = 1:M
        % Define training and test indices
        foldSize = floor(N / M);
        idx_test = (1+(k-1)*foldSize):(k*foldSize);
        idx_train = setdiff(1:N, idx_test);

        t_train = tspan(idx_train);
        t_test = tspan(idx_test);

        % Train model on training set
        initCond = [0; 0; 0.01; 0.01; 0.01];
        [t_train_sol, x_train_sol] = ode45(@(t,x) LyapSP(t,x,modelType), t_train, initCond);

        % Interpolate parameter estimates onto test set
        theta_hat_interp = interp1(t_train_sol, x_train_sol(:,3:5), t_test, 'linear', 'extrap');

        % Initialize test trajectories
        x1_test = zeros(length(t_test),1);
        x1_hat_test = zeros(length(t_test),1);
        x1_test(1) = 0; x1_hat_test(1) = 0;

        % Simulate dynamics using interpolated parameters
        for i = 2:length(t_test)
            dt = t_test(i) - t_test(i-1);
            x_prev = x1_test(i-1);
            xh_prev = x1_hat_test(i-1);
            t_curr = t_test(i-1);

            % Control input
            u = c * sin(d * t_curr) + f + 0.1 * sin(0.2*t_curr) + 0.5 * sin(0.3*t_curr);

            % True system update
            x1_test(i) = x_prev + dt * (-x_prev^3 + theta_true(1)*tanh(x_prev) + theta_true(2)/(1+x_prev^2) + u);

            % Extract interpolated parameters
            th = theta_hat_interp(i,:)';

            % Compute regression vector φ depending on the model
            switch modelType
                case 1
                    phi = [tanh(x_prev); 1/(1+x_prev^2); 0];
                case 2
                    phi = [x_prev^2; sin(x_prev); x_prev^2/(x_prev^2+1)];
                case 3
                    phi = [x_prev^3; tanh(x_prev); 1/(1+x_prev^2)];
                case 4
                    phi = [x_prev; x_prev^2; x_prev^3];
                case 5
                    phi = [exp(-x_prev^2); exp(-(x_prev-1)^2); exp(-(x_prev+1)^2)];
                otherwise
                    error('Invalid model type');
            end

            % Prediction update using interpolated parameters
            e = x_prev - xh_prev;
            x1_hat_test(i) = xh_prev + dt * (th(1)*phi(1) + th(2)*phi(2) + th(3)*phi(3) + u + thetam*e);
        end

        % Compute energy of validation error for this fold
        e_cv = x1_test - x1_hat_test;
        fold_energy(k) = trapz(t_test, e_cv.^2);
    end

    % Store average error across folds for this model
    cv_energy(modelType) = mean(fold_energy);
end

% Display cv results
fprintf('\n===== Cross-Validation Results =====\n');
for i = 1:5
    fprintf('Model %d: Avg CV Error = %.5f\n', i, cv_energy(i));
end
[~, best_cv] = min(cv_energy);
fprintf('\nBest model from CV: Model %d\n', best_cv);

%% ========== SYSTEM DYNAMICS WITH LYAPUNOV SERIES PARALLEL LAW ==========
function xxdot = LyapSP(t, xx, modelType)
    global c d f gamma thetam theta_true

    % States and estimates
    x1 = xx(1); x1_hat = xx(2);
    theta1 = xx(3); theta2 = xx(4); theta3 = xx(5);

    % Control input
    u = c * sin(d * t) + f + 0.1 * sin(0.2*t) + 0.5 * sin(0.3*t);
    
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
    xxdot(1) = -x1^3 + theta_true(1)*tanh(x1) + theta_true(2)/(1+x1^2) + u;          % True system
    xxdot(2) = theta1*phi(1) + theta2*phi(2) + theta3*phi(3) + u + thetam*e1;        % Estimated system
    xxdot(3) = gamma(1) * e1 * phi(1);                                               % θ1 estimate
    xxdot(4) = gamma(2) * e1 * phi(2);                                               % θ2 estimate
    xxdot(5) = gamma(3) * e1 * phi(3);                                               % θ3 estimate
end
