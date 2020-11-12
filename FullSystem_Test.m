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
ranges = 400:100:600;
iterations = 2;

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
        [~, ind] = min(abs(scenario.detection.detect_list.range - ranges(n)));
        range_out(n, m) = scenario.detection.detect_list.range(ind);
        
    end
end

%% Test Result Calculation and Visualization

% Mean and variance
range_mean = mean(range_out, 2, 'omitnan');
range_var = var(range_out, [], 2, 'omitnan');
range_error = range_mean - ranges;

% Plot results
figure;
plot(range_error);
hold on;
plot(range_var);
grid on;

%% Save and Package Resultant Data

% Run all end-of-process tasks
EndProcess












