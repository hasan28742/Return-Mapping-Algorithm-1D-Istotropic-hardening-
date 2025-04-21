clc
close all
clear all 

% Define Material Parameters
E = 100.0e3;      % Elastic Modulus (MPa)
sigma0 = 750;   % Initial Yield Strength (MPa)
H = 25.0e3;       % Plastic Modulus (MPa)
%H=0
Et = (E * H) / (E + H); % Tangent Modulus (MPa)

% Define Key Points from Calculation
%time_key = [1, 2, 4, 5];
time_key =[0,1 ,2];
%strain_key = [0, 0.02, -0.02, 0];
strain_key = [0, -.01,0]

% Plot strain history
figure(1)
plot(time_key, strain_key,'LineWidth',1.5)
hold on
title('provided stress-strain plot in the Problem');
grid on
xlabel('Time (arbitrary units)');
ylabel('Strain \epsilon (mm/mm)');

% Plot plastic hardeing modulus
figure(2)
example_plastic = linspace(0,0.1, 10)
example_stress = H*example_plastic + sigma0
plot(example_plastic, example_stress)
title('Plot plastic hardeing modulus');
xlabel('Plastic strain \epsilon (mm/mm)');
ylabel('Yield stress (MPa)');

% --- Discretize Loading Path ---
num_steps = 1000; % Total number of steps for the entire simulation
time = linspace(time_key(1), time_key(end), num_steps);
% Interpolate strain values at each time step
total_strain = interp1(time_key, strain_key, time);


% Initialization variables at zero
stress = zeros(1, num_steps);          % Stress history (GPa)
plastic_strain = zeros(1, num_steps);  % Plastic strain history
total_accum_plastic_strain = 0;              % Accumulated plastic strain (scalar)
yield_strength = sigma0;               % Current yield strength (scalar, starts at sigma0=500MPa)

tolerance = 1e-9; % Tolerance for yield check

%Go over each step, determine if yielding or not at each step
for i = 2:num_steps % Start from the second step

    % Current state variables (from previous step i-1)
    prev_stress = stress(i-1);
    prev_yield_strength = yield_strength; % Use the latest updated yield strength
    prev_accum_plastic_strain = total_accum_plastic_strain;
    prev_plastic_strain = plastic_strain(i-1);

    % Calculate total strain increment for this step
    delta_epsilon = total_strain(i) - total_strain(i-1);

    % Elastic stress redictor(assume elastic)
    sigma_trial = prev_stress + E * delta_epsilon;

    % Check Yield Condition (This is your yield surface)
    f_trial = abs(sigma_trial) - prev_yield_strength;

    % Make choice if elastic or plastic
    % Elastic- Yield strength and accumulated plastic strain don't change
    if f_trial <= tolerance
        % Elastic Step
        stress(i) = sigma_trial;
        plastic_strain(i) = prev_plastic_strain; % No change in plastic strain
    % Plastic- isotrpic hardening    
    else
         % Calculate plastic strain increment magnitude 
        delta_gamma = f_trial / (E + H);  %this is always positive since f_trial is positeive here

        % Calculate plastic strain increment, note it uses direction of
        % stress to add or subtract the plastic strain increment
        delta_eps_p = delta_gamma * sign(sigma_trial);

        % Update total plastic strain for this step
        plastic_strain(i) = prev_plastic_strain + delta_eps_p;

        % Update total accumulated plastic strain (magnitude for hardening)
        total_accum_plastic_strain = prev_accum_plastic_strain + delta_gamma;

        % Update new yield strength for the *next* step's check
        yield_strength = sigma0 + H * total_accum_plastic_strain

        % Update stress (correct back to the new yield surface)
        stress(i) = sigma_trial - E * delta_eps_p;
    end
end

% Stress History
figure(3);
plot(time, stress, '-b', 'LineWidth', 1.5);
xlabel('Time (arbitrary units)');
ylabel('Stress \sigma (MPa)');
title('stress vs strain plot')
grid on


% Stress-Strain Curve
figure(4)
plot(total_strain, stress, '-b', 'LineWidth', 1.5);
xlabel('Strain \epsilon');
ylabel('Stress \sigma (MPa)');
title('Stress-Strain Curve - Incremental Method 1.a');
grid on
 
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
 