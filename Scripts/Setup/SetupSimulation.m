%% MiRS - Simulation Setup
%{

    Sean Holloway
    MiRS Simulation Init File
    
    This file specifies simulation parameters for MiRS simulation.
    
%}

%% Simulation Parameter Setup

save_format.list = {'.png'};

% Radar simulation and processing setup
scenario.simsetup = struct( ...
    ...
    ... % Transceiver Trajectory Properties
    'radar_pos',    [0; 0; 0], ...              % Position of radar unit
    'radar_vel',    [0; 0; 0], ...              % Velocity of radar unit
    'radar_dir',    [0; 0], ...                 % Azimuth-Elevation normal of radar unit
    ...
    ... % Simulation Properties
    'num_frames',   1, ...                      % Number of radar frames to simulate
    'readout',      true, ...                   % Read out target data T/F
    ...
    'sim_rate',     2^0, ...                    % Rate to divide fast x slow time simulation
    'clear_cube',   false, ...
    'send_alert',   false, ...                  % Send email alert T/F
    'attach_zip',   false, ...
    'alert_address', 'sholloway@intellisenseinc.com', ...
    ...                                         % Email address for status updates
    'filename',     'Test_MiRS', ...            % Filename to save data as
    'save_format',  save_format, ...            % File types to save figures
    'save_figs',    false, ...                  % Save figures T/F
    'save_mat',     false, ...                  % Save mat file T/F
    'reduce_mat',   false);                     % Reduce mat file for saving




