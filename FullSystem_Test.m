%% MiRS - Full System w/ Multiple Runs for Testing
%{

    Sean Holloway
    MiRS (Millimeter Wave Ranging and Synchronization System)
    MATLAB Simulation & Processing

    This shell file runs successive scripts and gauges progress.
%}

%% Housekeeping, Timing, and Path Management

% Add current folders to path
addpath(genpath(pwd));

% Run all start-of-process tasks
StartProcess;

%% Loop Through Test Parameters

% Test parameters described here
ranges = 50:50:1000;
iterations = 20;

% Initialize outputs
range_out = zeros(length(ranges), iterations);
vel_out = range_out;
az_out = range_out;
snr_out = range_out;
snr_calc = zeros(length(ranges), 1);

for n = 1:length(ranges)
    
    for m = 1:iterations
        %% Initialize Scenario Object
        
        % Initialization
        scenario = RadarScenario;
        
        %% Setup Structures for Simulation
        
        % Set up simulation parameters
        SetupSimulation
        
        % Set up radar target
        SetupTarget
        
        % Set up transceiver and channel parameters
        SetupRadarScenario
        
        %% Modify Test Parameters
        
        % Test parameters input here
        scenario.target_list.pos = ranges(n) * [1; 0; 0];
        scenario.sim.target_plat.InitialPosition = ranges(n) * [1; 0; 0];
        
        %% Run Simulation & Signal Processing
        
        % Perform main processes of simulation, signal and data processing
        Main
        
        %% Collect Data From Test Loop
        
        % Output parameters here
        if scenario.detection.detect_list.num_detect > 0
            [~, ind] = min(abs(scenario.detection.detect_list.range - ranges(n)));
            range_out(n, m) = scenario.detection.detect_list.range(ind);
            vel_out(n, m) = scenario.detection.detect_list.vel(ind);
            az_out(n, m) = scenario.detection.detect_list.az(ind);
            snr_out(n, m) = scenario.detection.detect_list.SNR(ind);
        else
            range_out(n, m) = nan;
            vel_out(n, m) = nan;
            az_out(n, m) = nan;
            snr_out(n, m) = nan;
        end
        
        
    end
    
    % Save ideal SNR
    snr_calc(n) = CalculateSNR(scenario, scenario.target_list.rcs, sqrt(sum(scenario.target_list.pos.^2)));
    
end

%% Test Result Calculation and Visualization

% Save results
save('Data/MAT Files/ErrorTest.mat', 'range_out', 'vel_out', 'az_out', 'snr_out', 'snr_calc');

% Mean and variance
range_mean = mean(range_out, 2, 'omitnan');
range_var = var(range_out, [], 2, 'omitnan');
range_error = range_mean - ranges';

vel_var = var(vel_out, [], 2, 'omitnan');
vel_error = mean(vel_out, 2, 'omitnan');

az_var = var(az_out, [], 2, 'omitnan');
az_error = mean(az_out, 2, 'omitnan');

% Plot results
figure;
plot(range_error);
hold on;
plot(range_var);
grid on;

figure;
plot(vel_error);
hold on;
plot(vel_var);
grid on;

figure;
plot(az_error);
hold on;
plot(az_var);
grid on;


%% Save and Package Resultant Data

% Run all end-of-process tasks
EndProcess












