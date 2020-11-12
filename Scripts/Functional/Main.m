%% MiRS - Main Processing Loop
%{

    Sean Holloway
    MiRS Main Processing Loop
    
    This file specifies performs simulation, signal processing, detection,
    data processing, and results collection for MiRS system.
    
%}

%% Main Loop

% Start timing for estimation
timeStart(scenario);

for loop = 1:scenario.simsetup.num_frames
    
    scenario.flags.frame = loop;
    
    for cpi = 1:scenario.radarsetup.cpi_fr
        
        scenario.flags.cpi = cpi;
        
        %% Radar Simulation
        
        % Run simulation to retrieve fast time x slow time x rx-channel signal
        scenario = RadarSimulation(scenario);
        
        %% Signal Processing
        
        % Perform signal processing on received signal
        scenario.cube = SignalProcessing(scenario);
        
        % Generate 3-D Cartesian Mesh Grids
        %TODO: CHANGE THIS TO 3-COORDINATE SYSTEM
%         generateCoordinates(scenario);
        
        %% Single CPI Data Processing
        
        % Perform single-frame radar detection
        if scenario.simsetup.par_cfar
            scenario.detection = DetectionSingle_Parallel(scenario);
        else
            scenario.detection = DetectionSingle(scenario);
        end
        
        %% Loop Update Procedures
        
        % Read out CPI update
        CPIUpdate(scenario);
        
        % Read out estimated time of completion
        timeUpdate(scenario, 1, 'loops')
        
    end
    
    %% Multiple CPI Data Processing
    
    %TODO: UPDATE ALL
    %
    % Perform binary integration and coordinate determination
    scenario.detection = DetectionMultiple(scenario);
    
    % Read out detection data
    if scenario.simsetup.readout
        readOut(scenario);
    end
    
    % Save detection data
    saveMulti(scenario);
    
    % Update multi-target tracking system
%     scenario.multi = Tracking(scenario);
    %}
    
    %% Single Slice Visualization
    
    %{
    % View Range-Doppler heat map of center azimuth-elevation direction
    % viewRDCube(scenario, 'heatmap')
    
    % View Range-Doppler surface of center azimuth-elevation direction
    % viewRDCube(scenario, 'surface')
    
    % View Range-Angle heat map of zero-doppler slice
    % viewRACube(scenario, 'heatmap')
    
    % View Range-Angle surface of zero-doppler slice
    % viewRACube(scenario, 'surface')
    
    % View Range-Angle PPI of zero-doppler slice
    % viewRACube(scenario, 'PPI')
    
    % View single-frame detection heatmap
    % viewDetections(scenario, 'heatmap')
    
    % View single-frame detection PPI
    % viewDetections(scenario, 'PPI')
    %}
    
end

%% Data Visualization

% Display detections in 3D scatter plot
% viewDetections3D(scenario);

% Display result visualization plots
% viewTracking(scenario);








