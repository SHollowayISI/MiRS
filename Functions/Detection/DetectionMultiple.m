function [detection] = DetectionMultiple(scenario)
%DETECTIONMULTIPLE Performs multiple-CPI detection for MiRS project
%   Takes scenario object as input, provides scenario.detection object as
%   output, containing information about detected targets.

%% Unpack Variables

detection = scenario.detection;
radarsetup = scenario.radarsetup;
cube = scenario.cube;

%% Perform Binary Integration

% Find over-threshold detections
bw_cube = (detection.detect_cube_multi >= radarsetup.det_m);

% Average power for multiple-detection indices
avg_cube = bw_cube .* (detection.pow_cube_multi ./ detection.detect_cube_multi);
avg_cube(isnan(avg_cube)) = 0;

% Sum over angle information
rd_cube = sum(avg_cube, 3);

%% Determine Individual Object Coordinates

% Find connected objects in R-D cube
cc = bwconncomp(bw_cube);
regions = regionprops(cc, rd_cube, 'WeightedCentroid');

% Generate list of detection coordinates
detection.detect_list.range = [];
detection.detect_list.vel = [];
detection.detect_list.cart = [];
detection.detect_list.SNR = [];
detection.detect_list.num_detect = length(regions);

% Determine Centroid of azimuth-elevation slice
for n = 1:length(regions)
    
    % Store direct coordinates
    detection.detect_list.range(end+1) = interp1(cube.range_axis, regions(n).WeightedCentroid(2));
    detection.detect_list.vel(end+1) = interp1(cube.vel_axis, regions(n).WeightedCentroid(1));
    
    % Store derived coordinates
%     detection.detect_list.cart(:,end+1) = detection.detect_list.range(end) * ...
%         [cosd(detection.detect_list.el(end)) * cosd(detection.detect_list.az(end)); ...
%         cosd(detection.detect_list.el(end)) * sind(detection.detect_list.az(end));
%         sind(detection.detect_list.el(end))];
    
    % Store SNR
    detection.detect_list.SNR(end+1) = 10*log10(max(avg_cube, [], 'all')) ...
        - detection.noise_pow;
end



end

