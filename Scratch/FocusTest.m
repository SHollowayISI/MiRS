
%% Bookkeeping

close all;
clear variables;
tic;

%% Setup

% Signal properties
centerFrequencies = [0.5, 1, 3, 6] * 1e9;
bandwidth = 10e6;

% Coordinate setup
numberOfNodes = 10;
heightAboveTerrain = 1e3;
heightVariation = 10;
fieldSize = 1.6e3;
groundDistanceToBase = 10e3;

% Coordinate calculations
baseCoordinates = [groundDistanceToBase; 0; 0];
masterCoordinates = [0; 0; heightAboveTerrain];
slaveCoordinates = [fieldSize * (rand(2, numberOfNodes-1)-0.5);...
    heightAboveTerrain + heightVariation*(rand(1, numberOfNodes-1)-0.5)];
nodeCoordinates = [masterCoordinates, slaveCoordinates];

% Measurement calculations
positionVariance = 0;
measuredCoordinates = nodeCoordinates + positionVariance * randn(size(nodeCoordinates));

%% Focusing Algorithm

% Run focus algorithm
focusLoss = CalculateFocusLoss(baseCoordinates, nodeCoordinates, measuredCoordinates, centerFrequencies)
focusLossShifted = CalculateFocusLoss(baseCoordinates, nodeCoordinates, measuredCoordinates, centerFrequencies)

%% Function declarations

function [focusLoss, arrayGain] = CalculateFocusLoss(baseCoordinates, realCoordinates, measuredCoordinates, centerFrequencies)

% Calculate times of flight
[~, trueTimeDelay] = CalculateDelay(baseCoordinates, realCoordinates);
[~, measuredTimeDelay] = CalculateDelay(baseCoordinates, measuredCoordinates);

% Calculate timing error
timeError = (trueTimeDelay - measuredTimeDelay)';

% Calculate gain
phases = 2 * pi * centerFrequencies .* timeError;
amplitudes = sum(exp(1i * phases), 1);
arrayGain = 20*log10(abs(amplitudes));
focusLoss = (20*log10(size(realCoordinates, 2)) - arrayGain)';

end

function [timeOfFlight, timeDelay] = CalculateDelay(baseCoord, nodeCoords)

ranges = sqrt(sum((nodeCoords - baseCoord).^2, 1));

timeOfFlight = ranges / physconst('Lightspeed');
maxTOF = max(timeOfFlight);
timeDelay = timeOfFlight - maxTOF;

end


















