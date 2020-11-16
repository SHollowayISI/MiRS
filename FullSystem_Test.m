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
vel_in = 0;
bearing_in = 0;

range_var = 10;
vel_var = 0.1;
bearing_var = 10;

iterations = 50;

% Initialize outputs
range_out = zeros(length(ranges), iterations);
vel_out = range_out;
aoa_out = range_out;
snr_out = range_out;
calc_out = range_out;

for n = 1:length(ranges)
    
    for m = 1:iterations
        %% Initialize Scenario Object
        
        % Initialization
        scenario = RadarScenario;
        
        %% Modify Test Parameters
        
        % Add uncertainty to measurements
        current_range = ranges(n) + range_var * (rand(1)-0.5);
        current_vel = vel_in + vel_var * (rand(1)-0.5);
        current_bearing = bearing_in + bearing_var * (rand(1)-0.5);
        
        % Test parameters input here
        tgt_pos_in = current_range * [cosd(current_bearing); sind(current_bearing); 0];
        tgt_vel_in = current_vel * [cosd(current_bearing); sind(current_bearing); 0];
        
        %% Setup Structures for Simulation
        
        % Set up simulation parameters
        SetupSimulation
        
        % Set up radar target
        SetupTarget
        
        % Set up transceiver and channel parameters
        SetupRadarScenario
        
        %% Run Simulation & Signal Processing
        
        % Perform main processes of simulation, signal and data processing
        Main
        
        %% Collect Data From Test Loop
        
        % Output parameters here
        if scenario.detection.detect_list.num_detect > 0
            [~, ind] = min(abs(scenario.detection.detect_list.range - current_range));
            range_out(n, m) = scenario.detection.detect_list.range(ind) - current_range;
            vel_out(n, m) = scenario.detection.detect_list.vel(ind) - current_vel;
            aoa_out(n, m) = scenario.detection.detect_list.aoa(ind) - current_bearing;
            snr_out(n, m) = scenario.detection.detect_list.SNR(ind);
            calc_out(n, m) = CalculateSNR(scenario, scenario.target_list.rcs, sqrt(sum(scenario.target_list.pos.^2)));
        else
            range_out(n, m) = nan;
            vel_out(n, m) = nan;
            aoa_out(n, m) = nan;
            snr_out(n, m) = nan;
            calc_out(n, m) = nan;
        end
        
    end
   
end

%% Test Result Calculation and Visualization

% Save results
save('MAT Files/Data/Range_Error_Test_Final.mat', 'range_out', 'vel_out', 'aoa_out', 'snr_out', 'calc_out');


%% Save and Package Resultant Data

% Run all end-of-process tasks
EndProcess












