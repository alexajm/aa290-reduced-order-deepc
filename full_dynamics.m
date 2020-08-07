function [x_next, y] = full_dynamics(x, u)
% Load linearized matrices
global FOM process_cov measure_cov
A = FOM.A;
B = FOM.B;
C = FOM.C;

% Get next state, measurement
x_next = A*x + B*u + process_cov*randn(size(x));
y = C*x_next + measure_cov*randn(size(C,1),1);
end