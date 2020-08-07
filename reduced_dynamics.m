function [x_next, y] = reduced_dynamics(x, u)
% Load linearized matrices
global ROM
A = ROM.A;
B = ROM.B;
C = ROM.C;

% Get next state
x_next = A*x + B*u;
y = C*x_next;
end