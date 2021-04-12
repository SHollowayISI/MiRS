
clear variables;
close all;

% load('WithWeighting.mat');
% load('WithoutWeighting.mat');

filenames = {'WithWeighting', 'WithoutWeighting'};

field_ind = 1:11;
x_axis = 0:20:200;

leg = {};
for n = 1:length(x_axis)
    leg{n} = sprintf('Field Size %dm', x_axis(n));
end

for file = 1:length(filenames)
    
    load([filenames{file}, '.mat']);
    
    figure('Name', 'Inclination Error', 'units','normalized','outerposition',[0 0 1 1]);
    plot(x_axis, squeeze(rms(inc_out(:,field_ind,:), 1))');
    ylim([0 10]);
    grid on;
    xlabel('Master Node Offset', 'FontWeight', 'bold');
    ylabel('Inclination Error [deg]', 'FontWeight', 'bold');
    legend(leg)
    
    figure('Name', 'Yaw Error', 'units','normalized','outerposition',[0 0 1 1]);
    plot(x_axis, squeeze(rms(phi_out(:,field_ind,:), 1))');
    ylim([0 5]);
    grid on;
    xlabel('Master Node Offset', 'FontWeight', 'bold');
    ylabel('Yaw Error [deg]', 'FontWeight', 'bold');
    legend(leg)
    
    figure('Name', 'Total Position Error', 'units','normalized','outerposition',[0 0 1 1]);
    plot(x_axis, squeeze(mean(rms_out(:,field_ind,:), 1))');
    ylim([0 0.05]);
    grid on;
    xlabel('Master Node Offset', 'FontWeight', 'bold');
    ylabel('Total Position Error [m]', 'FontWeight', 'bold');
    legend(leg)
    
    figure('Name', 'Horizontal Position Error', 'units','normalized','outerposition',[0 0 1 1]);
    plot(x_axis, squeeze(mean(xy_rms_out(:,field_ind,:), 1))');
    ylim([0 0.005]);
    grid on;
    xlabel('Master Node Offset', 'FontWeight', 'bold');
    ylabel('Horizontal Position Error [m]', 'FontWeight', 'bold');
    legend(leg)
    
    figure('Name', 'Altitude Error', 'units','normalized','outerposition',[0 0 1 1]);
    plot(x_axis, squeeze(mean(z_rms_out(:,field_ind,:), 1))');
    ylim([0 0.05]);
    grid on;
    xlabel('Master Node Offset', 'FontWeight', 'bold');
    ylabel('Altitude Position Error [m]', 'FontWeight', 'bold');
    legend(leg)
    
    figure('Name', 'Timing Error', 'units','normalized','outerposition',[0 0 1 1]);
    plot(x_axis, (1e12 / physconst('Lightspeed')) * squeeze(mean(range_out(:,field_ind,:), 1))');
    ylim([0 100]);
    grid on;
    xlabel('Master Node Offset', 'FontWeight', 'bold');
    ylabel('Timing Error [ps]', 'FontWeight', 'bold');
    legend(leg)
    
    SaveFigures(filenames{file}, 'Figures/WeightingComparison', '.png');
    close all;
    
end


function [] = SaveFigures(save_name,fig_path, format)
%SAVEFIGURES Saves all open figures
%   Saves all open figures to fig_path directory, with name save_name plus
%   title of figure.

if ~exist(fig_path, 'dir')
    mkdir(fig_path)
end

FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    FigName   = get(FigHandle, 'Name');
    saveas(FigHandle, fullfile(fig_path, [save_name, ' ', FigName, format]));
end

% Display update to command window
disp([format, ' Figures saved in ', save_name]);

end





