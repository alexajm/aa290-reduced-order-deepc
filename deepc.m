% DeePC controller
%
% inputs
%       -x0: initial state
%       -traj: reference trajectory
%       -tf: time to simulate to
%       -dt: length of timestep
%       -H: struct of Hankels of offline data Up, Yp, Uf, Yf
%       -Np: backward time window (past)
%       -Nf: forward time window (future)
%
% outputs
%       -x: simulated state trajectory
%       -u: control policy
%
function [x,u] = deepc(x0, traj, tf, dt, H, Np, Nf)
% Globals
global FOM

% Setup
m = size(FOM.B,2);
n = size(FOM.A,2);
p = size(FOM.C,1);
lam_g = 1e2;
lam_y = 1e5;
L = Np + Nf;

% Load data
Up = H.Up;
Uf = H.Uf;
Yp = H.Yp;
Yf = H.Yf;

% Construct performance matrix
Qbase = 1e2 * ones(p,1);
Qdiag = repmat(Qbase,Nf,1);
Q = diag(Qdiag);

% Control constraints
% NOTE: currently hardcoded to large_synthetic
u_UB = .25;
u_LB = -.25;

% Run
x = x0;
u = [];
u_in = [];
y_out = FOM.C*x;
for t = 1:tf/dt
    % Get reference trajectory for this window
    ref = traj(p*(t-1)+1:p*(t+Nf-1));
   
    % CVX
    cvx_begin quiet
    cvx_precision default
    variables g(size(Up,2),1) ustar(Nf*m,1)
    variables y(Nf*p,1) sig(length(y_out),1)
    minimize (lam_g*norm(g,1) + lam_y*norm(sig,1) +...
              quad_form(y - ref, Q) +...
              quad_form(ustar, eye(length(ustar))))
    subject to
        [Up(1:length(u_in),:); Yp(1:length(y_out),:); Uf; Yf] * g ==...
            [u_in; y_out; ustar; y] +...
            [zeros(size(u_in)); sig; zeros(Nf*(m+p),1)];
        u_LB <= ustar; ustar <= u_UB;
        %ylim_LB <= y; y <= ylim_UB;
        sig >= 0;
    cvx_end
    ustar = Uf * g;
    
    % For when CVX is back on its Shit
    if isnan(ustar)
        ustar(1:m) = zeros(m,1);
    end
    
    % Update control sequence
    u = [u ustar(1:m)];
    
    % Progress dynamics
    [x_next, y_next] = full_dynamics(x(:,end),u(:,end));
    x = [x x_next];

    % Update input/output data
    u_in = [u_in(m+1:end); u(:,end)];
    y_out = [y_out(p+1:end); y_next];
end
end

