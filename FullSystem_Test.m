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
bearings = 0:5:90;
range_in = 250;
iterations = 10;

% Initialize outputs
range_out = zeros(length(bearings), iterations);
vel_out = range_out;
az_out = range_out;
aoa_out = range_out;
snr_out = range_out;
snr_calc = zeros(length(bearings), 1);

for n = 1:length(bearings)
    
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
        scenario.target_list.pos = range_in * [cosd(bearings(n)); sind(bearings(n)); 0];
        scenario.sim.target_plat.InitialPosition = range_in * [cosd(bearings(n)); sind(bearings(n)); 0];
        
        %% Run Simulation & Signal Processing
        
        % Perform main processes of simulation, signal and data processing
        Main
        
        %% Collect Data From Test Loop
        
        % Output parameters here
        if scenario.detection.detect_list.num_detect > 0
            [~, ind] = min(abs(scenario.detection.detect_list.range - range_in));
            range_out(n, m) = scenario.detection.detect_list.range(ind);
            vel_out(n, m) = scenario.detection.detect_list.vel(ind);
            az_out(n, m) = scenario.detection.detect_list.az(ind);
            aoa_out(n, m) = scenario.detection.detect_list.aoa(ind);
            snr_out(n, m) = scenario.detection.detect_list.SNR(ind);
        else
            range_out(n, m) = nan;
            vel_out(n, m) = nan;
            az_out(n, m) = nan;
            aoa_out(n, m) = nan;
            snr_out(n, m) = nan;
        end
        
        
    end
    
    % Save ideal SNR
    snr_calc(n) = CalculateSNR(scenario, scenario.target_list.rcs, sqrt(sum(scenario.target_list.pos.^2)));
    
end

%% Test Result Calculation and Visualization

% Save results
save('Data/MAT Files/AoA_Test.mat', 'range_out', 'vel_out', 'az_out', 'aoa_out', 'snr_out', 'snr_calc');

% Mean and variance
range_mean = mean(range_out, 2, 'omitnan');
range_var = var(range_out, [], 2, 'omitnan');
range_error = range_mean - range_in;

vel_var = var(vel_out, [], 2, 'omitnan');
vel_error = mean(vel_out, 2, 'omitnan');

az_mean = mean(az_out, 2, 'omitnan');
az_var = var(az_out, [], 2, 'omitnan');
az_error = az_mean - bearings';

aoa_mean = mean(aoa_out, 2, 'omitnan');
aoa_var = var(aoa_out, [], 2, 'omitnan');
aoa_error = aoa_mean - bearings';

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












