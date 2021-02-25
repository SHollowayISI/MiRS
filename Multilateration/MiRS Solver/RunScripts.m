
clear variables
close all;

tic;
field_sizes = [0.01, 0.1, 1, 10, 100];
offsets = [0, 1, 10, 100, 200];
iterations = 1000;
base_distances = [10e3, 30e3, 100e3];

theta_out = zeros(iterations, length(field_sizes), length(offsets));
psi_out = zeros(iterations, length(field_sizes), length(offsets));
phi_out = zeros(iterations, length(field_sizes), length(offsets));
ang_out = zeros(iterations, length(field_sizes), length(offsets));
rms_out = zeros(iterations, length(field_sizes), length(offsets));
xy_rms_out = zeros(iterations, length(field_sizes), length(offsets));
z_rms_out = zeros(iterations, length(field_sizes), length(offsets));
gains_out = zeros(4, length(base_distances), iterations, length(field_sizes), length(offsets));
focused_gains_out = zeros(4, length(base_distances), iterations, length(field_sizes), length(offsets));
focus_point_out = zeros(3, length(base_distances), iterations, length(field_sizes), length(offsets));

for f = 1:length(field_sizes)
    
    for o = 1:length(offsets)
        
        fprintf('Starting loop for field size %d and offset %d\n', field_sizes(f), offsets(o));
        toc;
        
        for i = 1:iterations
            
            if mod(i, 10) == 0
                fprintf('Starting Iteration: %d\n', i);
                toc;
            end
            
            field_size_in = field_sizes(f);
            base_distance_in = base_distances;
            offset_in = offsets(o);
            
            MiRS_Multilateration_3Dimension_Angle_v5;
            
            theta_out(i, f, o) = theta;
            psi_out(i, f, o) = psi;
            phi_out(i, f, o) = phi;
            ang_out(i, f, o) = total;
            rms_out(i, f, o) = err_rms;
            xy_rms_out(i, f, o) = err_xy_rms;
            z_rms_out(i, f, o) = err_z_rms;
            gains_out(:, :, i, f, o) = measuredGain;
            focused_gains_out(:, :, i, f, o) = gain;
            focus_point_out(:, :, i, f, o) = guessError;
            
        end
    end
end

%% Plotting

% close all;
%
% figure;
% c = [0.5, 1, 3, 6];
% for f = 1:3
%     for n = 1:4
%         scatter(1:iterations, gains_out(n,:,f),  '.');
%         hold on;
%     end
% end
% grid on;
% ylim([-20 20])
%
%
% figure;
% c = [0.5, 1, 3, 6];
% for f = 1:3
%     for n = 1:4
%         scatter(1:iterations, focused_gains_out(n,:,f),  '.');
%         hold on;
%     end
% end
% grid on;
% ylim([-20 20])

