function [scenario_out] = RadarSimulation(scenario)
%RADARSIMULATION Generates simulated radar response for MiRS
%   Takes scenario.sim, .simsetup, .radarsetup, .target_list as inputs and
%   outputs scenario.rx_sig struct containing the received signal.

%% Unpack Variables

scenario_out = scenario;
target_list = scenario.target_list;
sim = scenario.sim;
simsetup = scenario.simsetup;
radarsetup = scenario.radarsetup;
flags = scenario.flags;

%% Setup

% Physical constants
c = physconst('LightSpeed');

% Derived variables
num_ant = radarsetup.n_rx * radarsetup.n_tx;

%% Simulation

% Allocate size for signals
rx_sig = zeros(radarsetup.n_s - radarsetup.drop_s, radarsetup.n_p, num_ant);

% Generate the pulse
base_sig = sim.waveform();

% Transmit the pulse
tx_sig = sim.transmitter(base_sig);

for block = 1:radarsetup.n_p
    
    for chirp = 1:radarsetup.n_tx
        
        % Update target position and velocities
        [target_list.pos, target_list.vel] = sim.target_plat(radarsetup.t_ch);
        
        % Get range and angle to all targets
        [~, tgt_ang] = rangeangle(target_list.pos, simsetup.radar_pos);
        
        % Set up weights to implement TDM-MIMO
        weights = [0; 0];
        weights(chirp) = 1;
        
        % Radiate signal towards all targets
        sig = sim.radiator(tx_sig, tgt_ang, weights);
        
        % Propogate signal to the target through two-way channel
        sig = sim.channel(sig, ...
            simsetup.radar_pos, target_list.pos, ...
            simsetup.radar_vel, target_list.vel);
        
        % Reflect the pulse off of the target
        sig = sim.target(sig);
        
        % Collect the reflected target at the antenna array
        sig = sim.collector(sig, tgt_ang);
        
        % Apply receiver noise and gains
        sig = sim.receiver(sig);
        
        % Dechirp signal
        sig = dechirp(sig, base_sig);
        
        % Apply phase coding to result signal
        sig = reshape(sig, length(sig), 1, []);
        
        % Calculate indices for virtual antenna
        ch_ind = (1:radarsetup.n_rx) + radarsetup.n_rx * (chirp - 1);
        
        % Save Rx signal by fast time x slow time x virtual antenna
        rx_sig(:,block,ch_ind) = sig((radarsetup.drop_s+1):end,:,:);
        
    end
    
end


%% Re-pack Variables

scenario_out.rx_sig = rx_sig;

% Note: Do not export target list, since it is modified during


end