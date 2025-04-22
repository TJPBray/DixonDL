%% Shows different approaches to solve Ax = S when A and S are vectors and x is a scalar

% Define A 
A = [0 1 2 3 4 5 6 7]';

% Define ground truth x
x_gt = 10; 

% Define S (measured signal)

Snoisefree = x_gt*A;

% Add noise to lowest values 

noise = [10 0 0 0 0 0 0 0 ]';

Snoisy = Snoisefree + noise; 

%% Get x by element-wise division and then averaging
x0 = Snoisy./A;

x1 = mean(x0)

sse1 = sum ((x1*A - Snoisy).^2); 

%% Get x by pseudoinverse
x2 = pinv(A)*Snoisy

sse2 = sum ((x2*A - Snoisy).^2); 

