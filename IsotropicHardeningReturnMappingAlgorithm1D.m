clc
close all;

% Material properties
E = 100e3;            % Elastic modulus (MPa)
H = 30e3;             % Hardening modulus (MPa)
sigma0 = 500;         % Initial yield strength (MPa)

% Time steps from image
time_step = [1.0, 1.125, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0];


% Corresponding Strain history from handwritten table
strain_total = [0, 0.005, 0.01, 0.02, 0.01, 0.0, -0.01, -0.02, -0.01, 0]; 
N = length(strain_total);

% Initialize arrays
eps_e = zeros(1, N);  % Elastic strain
stress = zeros(1, N);          % final stress experienced in a step
sigma_trial = zeros(1, N);     % elastic trial stress
eps_p = zeros(1, N);           % plastic strain
eps_p_bar = zeros(1, N);       % accumulated plastic strain of effective plastic strain
yield_stress = zeros(1, N);    % current yield strength at a step

% Initial yield stress
yield_stress(1) = sigma0;

% Return Mapping Algorithm
for i = 2:N
    %calculate incrmental strain between two steps
    delta_eps = strain_total(i) - strain_total(i-1);

    % Elastic predictor incremental ( assume it is totally elastic step)
    sigma_trial(i) = stress(i-1) + E * delta_eps; %trial stress
    f_trial = abs(sigma_trial(i)) - yield_stress(i-1); % yield function

    if f_trial <= 1e-8  % f<=0
        % Elastic step
        stress(i) = sigma_trial(i); % this stess will be experienced by the metal < yield stress
        eps_p(i) = eps_p(i-1);  %no change in plastic strain
        eps_p_bar(i) = eps_p_bar(i-1); % no change in accumulated plastic strain
        yield_stress(i) = yield_stress(i-1); %no change in yield stress

    else %f>0
        % Plastic correction (return mapping)
        delta_gamma = f_trial / (E + H);
        delta_eps_p = delta_gamma * sign(sigma_trial(i));  %updated incremental plastic strain ( +ve or -ve both might be happened)

        eps_p(i) = eps_p(i-1) + delta_eps_p; %updated plastic strain unto this step ( +ve or -ve both might be happened)
        eps_p_bar(i) = eps_p_bar(i-1) + abs(delta_gamma);%updated accumulated plastic strain ( only +ve)
        yield_stress(i) = sigma0 + H * eps_p_bar(i); %new yield strenght due to accumulated plastic strain

        % Correct stress
        stress(i) = sigma_trial(i) - E * delta_eps_p; % stress experienced by the metal
    end

    %updated elastic strain after updating toatl and elastic strain
    eps_e(i) = strain_total(i) - eps_p(i); 

end

% Display results in table
T = table((time_step)', strain_total',  eps_e', eps_p', eps_p_bar', sigma_trial', stress', yield_stress',  ...
    'VariableNames', {'TimeStep', 'TotalStrain',  'ElasticStrain', 'PlasticStrain','AccumPlasticStrain'...
                      'TrialStress', 'UpdatedStress', 'YieldStress', });

disp(T);

% Plotting
figure;
subplot(2,1,1)
plot(time_step, stress, 'k-o', 'DisplayName', 'Updated Stress'); hold on;
plot(time_step, sigma_trial, 'c--x', 'DisplayName', 'Trial Stress');
plot(time_step, yield_stress, 'g-s', 'DisplayName', 'Yield Stress');
xlabel('Time Step'); ylabel('Stress (MPa)');
title('Stress and Yield Stress Over Steps');
legend('Location', 'best'); grid on;

subplot(2,1,2)
plot(time_step, strain_total, 'k-o', 'DisplayName', 'Total Strain'); hold on;
plot(time_step, eps_p, 'r-s', 'DisplayName', 'Plastic Strain');
plot(time_step, eps_p_bar, 'b-^', 'DisplayName', 'Accumulated Plastic Strain');
xlabel('Time_Step'); ylabel('Strain');
title('Strain Evolution Over Steps');
legend('Location', 'best'); grid on;
 

% 1. Time Step vs Stress
figure;
plot(time_step, stress, 'b-o', 'LineWidth', 1.5);
xlabel('Time Step');
ylabel('Stress (MPa)');
title('Time Step vs. Final Updated Stress');
grid on;

% 2. Time Step vs Strain Components
figure;
plot(time_step, strain_total, 'k-o', 'DisplayName', 'Total Strain'); hold on;
plot(time_step, eps_p, 'r-s', 'DisplayName', 'Plastic Strain');
plot(time_step, eps_e, 'm-^', 'DisplayName', 'Elastic Strain');
xlabel('Time Step');
ylabel('Strain');
title('Time Step vs. Strain Components');
legend('Location', 'best');
grid on;

% 3. Stress vs. Total Strain (Stress-Strain Curve)
figure;
plot(strain_total, stress, 'b-o', 'LineWidth', 2);
xlabel('Total Strain');
ylabel('Stress (MPa)');
title('Stress vs. Strain (total)');
grid on;
