%% MiRS - Example Radar Initialization File
%{

    Sean Holloway
    MiRS Init File
    
    This file specifies radar parameters for MiRS simulation.
    
%}

%% Radar Parameter Setup

% Radar simulation and processing setup
scenario.radarsetup = struct( ...
    ...
    ... % Waveform Properties
    'f_c',          78.5e9, ...             % Operating frequency in Hz
    'f_s',          40e6, ...               % ADC sample frequency in Hz
    't_ch',         1.0e-3, ...             % Chirp duration in seconds
    'bw',           5.0e9, ...              % Chirp bandwidth in Hz
    'n_p',          32, ...                  % Number of (MIMO) chirps per CPI
    'drop_s',       0, ...                  % Number of samples to drop
    'cpi_fr',       1, ...                  % Number of CPI per frame
    ...
    ... % Antenna Array Properties
    'n_tx',         2, ...                  % Number of elements in horizontal Tx array
    'd_tx',         2, ...                  % Distance between Tx elements in wavelengths
    'n_rx',         4, ...                  % Number of elements in horizontal Rx array
    'd_rx',         0.5, ...                % Distance between Rx elements in wavelengths
    ...
    ... % Transceiver Properties
    'tx_pow',       db2pow(13 - 30), ...    % Transmit power per antenna in Watts
    'rf_sys_loss',  0, ...                  % RF system loss in dB
    'rx_nf',        12, ...                 % Rx noise figure in dB
    ...
    ... % Antenna Properties
    'ant_gain',     db2pow(10), ...         % Antenna gain in absolute
    'ant_cos_pow',  [1 10], ...              % Cosine power of antenna pattern
    ...
    ... % Processing Properties
    'r_win',        'hanning', ...          % Window for range processing
    'd_win',        'blackmanharris', ...   % Window for doppler processing
    'az_win',       'hanning');             % Window for azimuth processing

%% Run Setup Scripts

% Set up Phased Array Toolbox system objects
scenario = PhasedSetup(scenario);





