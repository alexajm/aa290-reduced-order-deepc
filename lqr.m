function [x,u] = lqr(x0, traj, tf, dt)
% Globals
global FOM

% Setup
A = FOM.A;
B = FOM.B;
C = FOM.C;
m = size(FOM.B,2);
n = size(FOM.A,2);
p = size(FOM.C,1);
Q = 1e2 * ones(p,1);
R = 10 * ones(m,1);

% Get infinite-horizon Ricatti solution
[~,K,~] = idare(A,B,Q,R,[],[]);

u = [];
x = x0;
for t = 1:tf/dt
    % Progress dynamics
    u_next = -K*x(:,end);
    x_next = A*x(:,end) + B*u_next;
    
    % Save results
    u = [u u_next];
    x = [x x_next];
end
end