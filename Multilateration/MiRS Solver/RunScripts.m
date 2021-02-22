
tic;
% field_sizes = 0:10:200;
% field_sizes(1) = 1;
field_sizes = 0.1;
% offsets = 0:10:200;
offsets = 100;
iterations = 1000;

ang_out = zeros(iterations, length(field_sizes), length(offsets));
rms_out = zeros(iterations, length(field_sizes), length(offsets));
gains_out = zeros(iterations, 4, length(field_sizes), length(offsets));

for f = 1:length(field_sizes)
    for o = 1:length(offsets)
        
        fprintf('Starting loop for field size %d, offset %d\n', field_sizes(f), offsets(o));
        toc;
        
        for i = 1:iterations
            
            if mod(i, 10) == 0
                fprintf('Starting Iteration: %d\n', i);
                toc;
            end
            
            field_size_in = field_sizes(f);
            offset_in = offsets(o);
            
            MiRS_Multilateration_3Dimension_Angle_v3;
            
            ang_out(i, f, o) = total;
            rms_out(i, f, o) = err_rms;
            gains_out(i, :, f, o) = arrayGain;
            pos_gains_out(i,:,f,o) = rotatedArrayGain;
            
        end
    end
end

%% Plotting

close all;

figure;
c = [0.5, 1, 3, 6];
for n = 1:4
    scatter(1:iterations, gains_out(:,n),  '.');
    hold on;
end
grid on;
ylim([-20 20])


figure;
c = [0.5, 1, 3, 6];
for n = 1:4
    scatter(1:iterations, pos_gains_out(:,n),  '.');
    hold on;
end
grid on;
ylim([-20 20])

