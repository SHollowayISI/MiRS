
%% Bookkeeping

% clear variables;
% close all;
% tic;

%% User Inputs

% Simulation variables
num_nodes = 10;                         % Total number of nodes, including master
field_size_xy = 1600;                   % Width of total field of nodes, in meters
field_size_z = field_size_in;                       % Maximum separation in altitude between nodes
offset_z = offset_in;                         % Height of master node above slave nodes
error_var_r = range_error_in;                   % Standard deviation of range measurements
error_var_a = 1;                        % Standard deviation of angle measurements
error_var_z = min(0.1, field_size_z);   % Standard deviation of altitude measurements
digitalPhase = true;                    % Discretize phase of pointing algorithm

% Solver variables
max_iter = 10;                          % Maximum iterations to use for solver
min_iter = 1;                           % Minimum iterations to use for solver
r_factor = 1e0;                         % Extra weighting factor for range
a_factor = 1e0;                         % Extra weighting factor for angle
z_factor = 1e9;                         % Extra weighting factor for altitude
r_weight = r_factor/(error_var_r)^4;    % Covariance matrix of range measurements
a_weight = a_factor/(error_var_a)^2;    % Covariance matrix of angle measurements
z_weight = min(z_factor/(error_var_z)^2, 1e20);     %matrix of altitude measurements
cost_thresh = -Inf;                     % Log-cost function value required to break loop (unused)
cost_diff_thresh = -10;                 % Log-cost difference required to break loop (unused)
lambda_init = 0.01;                     % Initial lambda value for levenberg-marquadt
L_up = 11;                              % Lambda adjustment factor for LM algorithm
L_down = 9;                             % Lambda adjustment factor for LM algorithm
epsilon = 0.9;                          % Lambda adjustment factor for LM algorithm

%% Setup

% Calculations
num_var = 3*(num_nodes - 1);
num_eq = (num_nodes - 1) * num_nodes / 2;

% Initialize loop variables
cost = zeros(max_iter, 1);
lambda = zeros(max_iter, 1);
x_guess = repmat({zeros(num_var, 1)}, max_iter, 1);

% Indices
n = floor((3 + sqrt(8*(1:num_eq)' - 7))/2);
m = (1:num_eq)' - (n - 2) .* (n - 1)/2;

% Simulation Setup
x_real = field_size_xy * (rand(num_var, 1) - 0.5);
x_real(mod(1:length(x_real), 3) == 0) = field_size_z * (rand((num_nodes-1),1) - 0.5);
x_real = [0; 0; offset_z; x_real];

% Calculate exact ranges (squared)
n_ind = 3*floor((n-1)) + (1:3) - 3;
m_ind = 3*floor((m-1)) + (1:3) - 3;

% Calculate exact angles and ranges
r_sq = sum((x_real(n_ind + 3) - x_real(m_ind + 3)).^2, 2);
ang = atan2d(x_real(n_ind(:,2) + 3) - x_real(m_ind(:,2) + 3), x_real(n_ind(:,1) + 3) - x_real(m_ind(:,1) + 3));

% First guess using blurry real points from master node
guess_ang = ang(m == 1) + error_var_a * randn(size(ang(m == 1)));
guess_r = sqrt(r_sq(m == 1)) + error_var_r * randn(size(r_sq(m == 1)));
guess_x = guess_r .* cosd(guess_ang);
guess_y = guess_r .* sind(guess_ang);
z_real = x_real(mod(1:length(x_real),3) == 0);
z_real(1) = [];
guess_z = z_real + error_var_z * randn(length(z_real), 1);
x_guess{1} = reshape(([guess_x, guess_y, guess_z])', [], 1);

% Introduce error into range and angle measurement
r_sq = (sqrt(r_sq) + error_var_r * randn(size(r_sq))).^2;
ang = ang + error_var_a*randn(size(ang));

% Vectorize data
y = [r_sq; ang; guess_z];

%% Initial step

% Matrix setup
W = diag([r_weight*ones(num_eq, 1); a_weight*ones(num_eq, 1); z_weight*ones(num_nodes-1, 1)]);
J = jacobian(x_guess{1}, n, m, n_ind, m_ind, num_var, num_eq, offset_z);

% Evaluate cost function
y_hat = angleRangeFunctions(x_guess{1}, n, m, n_ind, m_ind, offset_z, z_real);
cost(1) = costFunction(y_hat, y, W, J, zeros(num_var, 1));

% Initialize lambda parameter
lambda(1) = lambda_init;

%% Main loop

for iter = 2:max_iter
    
    % Get Jacobian matrix
    J = jacobian(x_guess{iter-1}, n, m, n_ind, m_ind, num_var, num_eq, offset_z);
    
    % Levenberg-Marquadt Step
    h = LMStep(y_hat, y, W, J, lambda(iter-1));
    
    % Lambda iteration
    if LMMetric(y_hat, y, W, J, h, lambda(iter-1)) > epsilon
        cost(iter) = costFunction(y_hat, y, W, J, h);
        x_guess{iter} = x_guess{iter-1} + h;
        lambda(iter) = max(lambda(iter-1)/L_down, 1e-7);
    else
        cost(iter) = cost(iter-1);
        x_guess{iter} = x_guess{iter-1};
        lambda(iter) = min(lambda(iter-1)*L_up, 1e7);
    end
    
    % Update
    y_hat = angleRangeFunctions(x_guess{iter}, n, m, n_ind, m_ind, offset_z, z_real);
    
end

%% Error estimation & Plotting

% Degree of freedom reduction and error estimation
x_result = reshape([0; 0; offset_z; x_guess{end}], 3, num_nodes)';
x_exact = reshape(x_real, 3, num_nodes)';

% Kabsch algorithm for rotation
[U, r, lrms] = Kabsch(x_result', x_exact');
x_rotated = (U * x_result' + r)';

% Calculate rotation angles
theta = -asind(U(3,1));
psi = atan2d(U(3,2)/cosd(theta), U(3,3)/cosd(theta));
phi = atan2d(U(2,1)/cosd(theta), U(1,1)/cosd(theta));
total = acosd((trace(U)-1)/2);

% Calculate error
error_list = x_rotated - x_exact;
err_off = sqrt(sum(mean(error_list,1).^2));
err_rms = rms(sqrt(sum((error_list - mean(error_list,1)).^2, 2)));
err_xy_off = sqrt(sum(mean(error_list(:,1:2),1).^2));
err_xy_rms = rms(sqrt(sum((error_list(:,1:2) - mean(error_list(:,1:2),1)).^2, 2)));
err_z_off = mean(error_list(:,3),1);
err_z_rms = rms(error_list(:,3) - mean(error_list(:,3),1));

range_error = mean(sqrt(sum((error_list(1,:) - error_list(2:num_nodes,:)).^2, 2)));

% Info printout
%{
fprintf('\nSimulation Complete.\n');
fprintf('Angle Error: (%0.3f x %0.3f x %0.3f = %0.3f) [deg]\n', theta, phi, psi, total);
fprintf('\nTotal Distance: \nError Offset: %d [m], \nError RMS: %d [m]\n', err_off, err_rms);
fprintf('\nXY-Plane: \nError Offset: %d [m], \nError RMS: %d [m]\n', err_xy_off, err_xy_rms);
fprintf('\nAltitude: \nError Offset: %d [m], \nError RMS: %d [m]\n', err_z_off, err_z_rms);

% Plotting
close all;

% Plot points
figure;
scatter3(x_result(:,1), x_result(:,2), x_result(:,3), '.', 'r');
hold on;
scatter3(x_exact(:,1), x_exact(:,2), x_exact(:,3));
hold on;
scatter3(x_rotated(:,1), x_rotated(:,2), x_rotated(:,3), '.', 'b');
% hold on;
% scatter3(x_rotated_2(:,1), x_rotated_2(:,2), x_rotated_2(:,3), 'g');
grid on;
xlim([-1 1] * field_size_xy/2)
xlabel('x')
ylim([-1 1] * field_size_xy/2)
ylabel('y')
% zlim([-field_size_z field_size_z + offset_z])
zlim([-1 1] * field_size_xy/2)
%}

%% Focusing algorithm

baseDistance = [10e3, 30e3, 100e3];
% baseDistance = [10e3];
delta = 1e-12;

for n = 1:length(baseDistance)
    % Coherent gain calculation
    baseCoord = [baseDistance(n); 0; -1e3];
    centerFrequencies = [0.5, 1, 3, 6] * 1e9;
    
    % Focus point estimation
    rotatedBase = U \ (baseCoord - r);
    [~, measuredToF] = CalculateDelay(baseCoord, x_result');
    [~, rotatedToF] = CalculateDelay(rotatedBase, x_result');
    [~, exactToF] = CalculateDelay(rotatedBase, x_exact');
    [~, trueToF] = CalculateDelay(baseCoord, x_exact');
    
    deltaT = (trueToF - measuredToF);
    if digitalPhase
        deltaT = delta * round(deltaT / delta);
    end
%     for i = 1:length(centerFrequencies)
%         measuredGain(i,n) = pulseGain(deltaT', centerFrequencies(i));
%     end
    phases = 2 * pi * centerFrequencies' .* deltaT;
    amplitude = sum(exp(1i * phases), 2);
    measuredGain(:,n) = 20*log10(abs(amplitude)); 
    
    deltaT = (trueToF - rotatedToF);
    if digitalPhase
        deltaT = delta * round(deltaT / delta);
    end
%     for i = 1:length(centerFrequencies)
%         gain(i,n) = pulseGain(deltaT', centerFrequencies(i));
%     end
    phases = 2 * pi * centerFrequencies' .* deltaT;
    amplitude = sum(exp(1i * phases), 2);
    gain(:,n) = 20*log10(abs(amplitude));
    
    guessError(:,n) = rotatedBase - baseCoord;
    guessPhi = atan2d(guessError(2), baseCoord(1));
    guessTheta = atan2d(guessError(3), baseCoord(1));
end

%% Draw heatmap of gain
%{
freq = 1e9;
baseDistance = 10e3;
hat = 1e3;
maxDiff = 1000;
step = 10;
digitalPhase = true;
delta = 1e-12;

xVal = (baseDistance - maxDiff):step:(baseDistance + maxDiff);
yVal = (-maxDiff):step:(maxDiff);

gainMap = zeros(length(xVal), length(yVal));

for x_n = 1:length(xVal)
    for y_n = 1:length(yVal)
        
        evalPoint = [xVal(x_n); yVal(y_n); -hat];
        
        [~, measuredToF] = CalculateDelay(evalPoint, x_result');
        [~, trueToF] = CalculateDelay(evalPoint, x_exact');
        
        deltaT = (trueToF - measuredToF);
        if digitalPhase
            deltaT = delta * round(deltaT / delta);
        end
        gainMap(x_n, y_n) = pulseGain(deltaT', freq) - 10;
%         phases = 2 * pi * freq .* deltaT;
%         amplitude = sum(exp(1i * phases), 2);
%         gainMap(x_n, y_n) = 20*log10(abs(amplitude)) - 10;
        
    end
end

close all;
figure;
imagesc(yVal, xVal, gainMap')
set(gca, 'YDir', 'normal')
caxis([-10 10])
colorbar
hold on;
scatter(0, baseDistance, '.', 'r')
%}


%% Focusing Algorithm Functions

function gain = ...
    CalculateFocusLoss(baseCoordinates, realCoordinates, measuredCoordinates, centerFrequencies)

% Calculate times of flight
trueTimeOfFlight = CalculateDelay(baseCoordinates, realCoordinates);
measuredTimeOfFlight = CalculateDelay(baseCoordinates, measuredCoordinates);

% Calculate gain
phases = 2 * pi * centerFrequencies .* (trueTimeOfFlight - measuredTimeOfFlight)';
amplitudes = sum(exp(1i * phases), 1);
gain = 20*log10(abs(amplitudes))';

end

function [timeOfFlight, timeDelay] = CalculateDelay(baseCoord, nodeCoords)

ranges = sqrt(sum((nodeCoords - baseCoord).^2, 1));

timeOfFlight = ranges / physconst('Lightspeed');
maxTOF = max(timeOfFlight);
timeDelay = timeOfFlight - maxTOF;

end

function gain = pulseGain(timeDifferences, frequencies)

tau = timeDifferences - mean(timeDifferences);
f = frequencies;
t = (0:1:30000) * 1e-12;

x = exp(1i * (t - tau) * 2 * pi * f) .*  rect(10000 * 1e-12, 20000 * 1e-12, t-tau);
sumSig = sum(x, 1)';
gain = 10*log10(3*sum(abs(sumSig).^2) / length(sumSig));

end

function y = rect(a, b, x)

y = zeros(size(x));
y((x < b) & (x > a)) = 1;
y((x == a) | (x == b)) = 0.5;

end

%% Cost function
function y_hat = angleRangeFunctions(x_in, n, m, n_ind, m_ind, offset_z, z_real)

% Evaluate cost function
r_sq_est = zeros(length(n),1);
for i = 1:length(n)
    if (m(i) == 1)
        r_sq_est(i) = sum(x_in(n_ind(i,1:2)).^2) + (x_in(n_ind(i,3)) - offset_z).^2;
    else
        r_sq_est(i) = sum((x_in(n_ind(i,:)) - x_in(m_ind(i,:))).^2);
    end
end

% Cost function for angle
ang_est = zeros(length(n),1);
for i = 1:length(n)
    if (m(i) == 1)
        ang_est(i) = atan2d(x_in(n_ind(i,2)), x_in(n_ind(i,1)));
    else
        ang_est(i) = atan2d(x_in(n_ind(i,2)) - x_in(m_ind(i,2)), x_in(n_ind(i,1)) - x_in(m_ind(i,1)));
    end
end

% Cost function for z range
z_est = zeros(length(x_in)/3, 1);
for row = 1:length(z_est)
    z_est(row) = x_in(3*row) - z_real(row);
end

y_hat = [r_sq_est; ang_est; z_est];

end

%% Jacobian
function J = jacobian(x_in, n, m, n_ind, m_ind, num_var, num_eq, offset_z)

% Calculate squared range between points
r_sq_est = zeros(length(n),1);
for i = 1:length(n)
    if (m(i) == 1)
        r_sq_est(i) = sum(x_in(n_ind(i,1:2)).^2) + (x_in(n_ind(i,3)) - offset_z).^2;
    else
        r_sq_est(i) = sum((x_in(n_ind(i,:)) - x_in(m_ind(i,:))).^2);
    end
end

% Calculate jacobian
J = zeros(2*num_eq  + num_var/3, num_var);
for row = 1:num_eq
    if m(row) == 1
        J(row, n_ind(row, 1:2)) = 2 * (x_in(n_ind(row,1:2)));
        J(row, n_ind(row, 3)) = 2 * (x_in(n_ind(row,3)) - offset_z);
        J(row + num_eq, n_ind(row,1)) = -x_in(n_ind(row,2));
        J(row + num_eq, n_ind(row,2)) = x_in(n_ind(row,1));
    else
        J(row, n_ind(row,:)) = 2 * (x_in(n_ind(row,:)) - x_in(m_ind(row,:)));
        J(row, m_ind(row,:)) = 2 * (x_in(m_ind(row,:)) - x_in(n_ind(row,:)));
        J(row + num_eq, n_ind(row,1)) = x_in(m_ind(row,2)) - x_in(n_ind(row,2));
        J(row + num_eq, m_ind(row,1)) = x_in(n_ind(row,2)) - x_in(m_ind(row,2));
        J(row + num_eq, n_ind(row,2)) = x_in(n_ind(row,1)) - x_in(m_ind(row,1));
        J(row + num_eq, m_ind(row,2)) = x_in(m_ind(row,1)) - x_in(n_ind(row,1));
    end
end
J((num_eq+1):(2*num_eq),:) = rad2deg(J((num_eq+1):(2*num_eq),:) ./ r_sq_est);

% Calculate jacobian for z-constraint
for row = 1:(num_var/3)
    J(2*num_eq + row, 3*row) = 1;
end

end

%% Cost function
function chi_sq = costFunction(y_hat, y, W, J, h)

% Calculate chi-squared cost function value
chi_sq = (y' * W * y) + (y_hat' * W * y_hat) - 2 * (y' * W * y_hat) + ...
    (h' * J' * W * J * h) - 2 * ((y - y_hat)' * W * J * h);

end

%% Levenberg-Marquadt Step
function h = LMStep(y_hat, y, W, J, lambda)

% Calculate step
h = ((J' * W * J) + lambda * diag(J' * W * J) .* eye(size(J, 2))) \ (J' * W * (y - y_hat));


end

%% Lambda Update Metric
function rho = LMMetric(y_hat, y, W, J, h, lambda)

% Calculate metric for L-M parameter update
rho = (costFunction(y_hat, y, W, J, zeros(size(h))) - costFunction(y_hat, y, W, J, h)) ./ ...
    (h' * (lambda * diag(J' * W * J) .* eye(size(J, 2)) * h + (J' * W * (y - y_hat))));

end

%% Kabsch Algorithm
function[U, r, lrms] = Kabsch(P, Q, m)
	sz1 = size(P) ;
	sz2 = size(Q) ;
	if (length(sz1) ~= 2 || length(sz2) ~= 2)
		error 'P and Q must be matrices' ;
	end
	if (any(sz1 ~= sz2))
		error 'P and Q must be of same size' ;
	end
	D = sz1(1) ;         % dimension of space
	N = sz1(2) ;         % number of points
	if (nargin >= 3)
		if (~isvector(m) || any(size(m) ~= [1 N]))
			error 'm must be a row vector of length N' ;
		end 
		if (any(m < 0))
			error 'm must have non-negative entries' ;
		end
		msum = sum(m) ;
		if (msum == 0)
			error 'm must contain some positive entry' ;
		end
		m = m / msum ;     % normalize so that weights sum to 1
	else                 % m not supplied - use default
		m = ones(1,N)/N ;
	end
	p0 = P*m' ;          % the centroid of P
	q0 = Q*m' ;          % the centroid of Q
	v1 = ones(1,N) ;     % row vector of N ones
	P = P - p0*v1 ;      % translating P to center the origin
	Q = Q - q0*v1 ;      % translating Q to center the origin
	% C is a covariance matrix of the coordinates
	% C = P*diag(m)*Q' 
	% but this is inefficient, involving an N*N matrix, while typically D << N.
	% so we use another way to compute Pdm = P*diag(m),
	% which is equivalent to, but more efficient than,
	% Pdm = zeros(D,N) ;
	% for i=1:N
	% 	Pdm(:,i) = m(i)*P(:,i) ;
	% end
	Pdm = bsxfun(@times,m,P) ;
	C = Pdm*Q' ; 	
%	C = P*Q' / N ;       % (for the non-weighted case)       
	[V,S,W] = svd(C) ;   % singular value decomposition
	I = eye(D) ;
	if (det(V*W') < 0)   % more numerically stable than using (det(C) < 0)
		I(D,D) = -1 ;
	end
	U = W*I*V' ;
	r = q0 - U*p0 ;
	Diff = U*P - Q ;     % P, Q already centered
%	lrms = sqrt(sum(sum(Diff.*Diff))/N) ; % (for the non-weighted case)
	% to compute the lrms, we employ an efficient method, equivalent to: 
	% lrms = 0 ;
	% for i=1:N
	% 	lrms = lrms + m(i)*Diff(:,i)'*Diff(:,i) ;
	% end
	% lrms = sqrt(lrms) ;
	lrms = sqrt(sum(sum(bsxfun(@times,m,Diff).*Diff))) ; 
end

%% Error Functions

% function sigma = errorSigma(r_sq, 







