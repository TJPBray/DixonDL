%% Shows different approaches to solve Ax = S when A and S are vectors and x is a scalar

%% Extends v1 to do multiple voxels simultaneously


% Define A (common across voxels)
A = [0 1 2 3 4 5 6 7]';

% Define ground truth x (now a column vector)
x_gt = [10 20 30 70]; 

% Define S (measured signal) - now a matrix

Snoisefree = A*x_gt;

% Add noise to lowest values 

noise = randn(size(Snoisefree));

Snoisy = Snoisefree + noise; 


%% Get x by pseudoinverse
x_hat = pinv(A)*Snoisy

sse2 = sum ((x2*A - Snoisy).^2); 

