%% MiRS - Full System
%{

    Sean Holloway
    MiRS (Millimeter Wave Ranging and Synchronization System)
    MATLAB Simulation & Processing

    This shell file runs successive scripts and gauges progress.

    TODO:
    - SetupRadarScenario:
        + Update fields and values for new processing chain:
            * n_p:                      Number of chirps per CPI
            * drop_s:                   Number of samples to drop
            * r_win, d_win, az_win:     Window types
    - Main:
        + Update all files from PAN-UAS version
        + GenerateCoordinates:
            * Remove or modify this
        + DetectionSingle:
            * CFAR/Threshold/other?
            * Run CFAR for every azimuth swath?
        + DetectionMultiple:
            * New detection method for MiRS
            * Remove elevation dimension
        + Visualization methods:
            * Remove elevation swaths
    - RadarScenario:
        + Go through each method and remove/modify

    NEXT STEPS:
    - Modify to multi-static version
        + Associate targets with transceivers
        + Implement target trajectory
%}

%% Housekeeping, Timing, and Path Management

% Add current folders to path
addpath(genpath(pwd));

% Run all start-of-process tasks
StartProcess;

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

%% Run Simulation & Signal Processing

% Perform main processes of simulation, signal and data processing
Main

%% Save and Package Resultant Data

% Run all end-of-process tasks
EndProcess














