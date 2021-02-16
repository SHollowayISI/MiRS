
tic;
% field_sizes = 0:10:200;
% field_sizes(1) = 1;
field_sizes = 10;
% offsets = 0:10:200;
offsets = 100;
iterations = 10000;

ang_out = zeros(iterations, length(field_sizes), length(offsets));
rms_out = zeros(iterations, length(field_sizes), length(offsets));

for f = 1:length(field_sizes)
    for o = 1:length(offsets)
        
        fprintf('Starting loop for field size %d, offset %d\n', field_sizes(f), offsets(o));
        toc;
        
        for i = 1:iterations
            
            field_size_in = field_sizes(f);
            offset_in = offsets(o);
            
            MiRS_Multilateration_3Dimension_Angle_v3;
            
            ang_out(i, f, o) = total;
            rms_out(i, f, o) = err_rms;
            
        end
    end
end
