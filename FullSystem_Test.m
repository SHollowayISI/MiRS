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
vels = 0.05:0.05:0.2;
range_in = 250;
iterations = 10;

% Initialize outputs
range_out = zeros(length(vels), iterations);
vel_out = range_out;
az_out = range_out;
aoa_out = range_out;
snr_out = range_out;
snr_calc = zeros(length(vels), 1);

for n = 1:length(vels)
    
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
        scenario.target_list.pos = range_in * [1; 0; 0];
        scenario.sim.target_plat.InitialPosition = range_in * [1; 0; 0];
        scenario.target_list.vel = vels(n) * [1; 0; 0];
        scenario.sim.target_plat.Velocity = vels(n) * [1; 0; 0];
        
        %% Run Simulation & Signal Processing
        
        % Perform main processes of simulation, signal and data processing
        Main
        
        %% Collect Data From Test Loop
        
        % Output parameters here
        if scenario.detection.detect_list.num_detect > 0
            [~, ind] = min(abs(scenario.detection.detect_list.range - range_in));
            range_out(n, m) = scenario.detection.detect_list.range(ind);
            vel_out(n, m) = scenario.detection.detect_list.vel(ind);
            aoa_out(n, m) = scenario.detection.detect_list.aoa(ind);
            snr_out(n, m) = scenario.detection.detect_list.SNR(ind);
        else
            range_out(n, m) = nan;
            vel_out(n, m) = nan;
            aoa_out(n, m) = nan;
            snr_out(n, m) = nan;
        end
        
        
    end
    
    % Save ideal SNR
    snr_calc(n) = CalculateSNR(scenario, scenario.target_list.rcs, sqrt(sum(scenario.target_list.pos.^2)));
    
end

%% Test Result Calculation and Visualization

% Save results
<<<<<<< Updated upstream
save('MAT Files/Data/AoA_Test.mat', 'range_out', 'vel_out', 'az_out', 'aoa_out', 'snr_out', 'snr_calc');
=======
save('MAT Files/Data/Vel_vs_AoA_Test.mat', 'range_out', 'vel_out', 'az_out', 'aoa_out', 'snr_out', 'snr_calc');
>>>>>>> Stashed changes

% Mean and variance
aoa_mean = mean(aoa_out, 2, 'omitnan');
aoa_var = var(aoa_out, [], 2, 'omitnan');


%% Save and Package Resultant Data

% Run all end-of-process tasks
EndProcess












