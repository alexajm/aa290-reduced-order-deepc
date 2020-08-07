function x = koopman(x0, traj, tf, dt)
% Globals
global ROM FOM u_LB u_UB

% Setup
lam = 1;
Td = 1e4;
Ts = 1;
Tu = 5*Ts;
ulim = 100;
m = size(FOM.B,2);
n = size(FOM.A,2);
p = size(FOM.C,1);
sig = 0;
Nh = 10;

% Collect offline data
x = zeros(n,1);
a = [];
u = zeros(m,1);
u_dat = u;
for t = 1:Td
    % Change control
    if mod(t,Tu) == 0, u = 2*ulim*rand(m,1) - ulim; end
    
    % Apply control  
    [x_next, y_next] = full_dynamics(x, u);
    x = x_next;
    
    % Collect sample
    if mod(t,Ts) == 0
        a = [a y_next];
        u_dat = [u_dat u];
    end
end
b = a(:,m+1:end);
a = a(:,1:end-m);

% Lift data (monomials up to 4th order)
K = size(a,2);
psi_a = [];
psi_b = [];
%lift = @(m,k) [m(:,k)' m(:,k)'.^2 m(:,k)'.^3 m(:,k)'.^4]';
lift = @(m,k) [m(:,k)']';
for k = 1:K
    psi_a = [psi_a; lift(a,k)'];
    psi_b = [psi_b; lift(b,k)'];
end
N = size(psi_a,2);

% Combine lifted data, inputs
alpha = @(k) [lift(a,k); u_dat(:,k)];
beta = @(k) [lift(b,k); u_dat(:,k)];
gam_a = [];
gam_b = [];
for k = 1:K
    gam_a = [gam_a; alpha(k)'];
    gam_b = [gam_b; beta(k)'];
end

% Reshape data and inputs
gam_b_vec = reshape(gam_b, K*(N+m), 1);
gam_a_blk = [];
for row = 1:size(gam_a,2)
    gam_a_blk = blkdiag(gam_a_blk,gam_a);
end

% Approximate Koopman operator w/ lasso
U = lasso(gam_a_blk,gam_b_vec,'Lambda',lam);
U = reshape(U,N+m,N+m);

% Extract model matrices
A = U(1:N,1:N)';
B = U(N+1:end,1:N)';

% Identify projection operator
omega_func = @(k) A*lift(a,k) + B*u_dat(:,k);
omega = [];
for k = 1:K
    omega = [omega; omega_func(k)'];
end
P = (pinv(omega) * psi_b)';

% Get projected model matrices
Ahat = P*A;
Bhat = P*B;
Ahat_blk = [];
Bhat_blk = [];
for i = 1:Nh
    Ahat_blk = blkdiag(Ahat_blk,Ahat);
    Bhat_blk = blkdiag(Bhat_blk,Bhat);
end

% Koopman MPC
x = x0;
y = FOM.C * x;
u = [];
G = eye((Nh+1)*N);
g = ones((Nh+1)*N,1);
H = eye(Nh*m);
h = ones(Nh*m,1);
ustar_LB = repmat(u_LB,Nh,1);
ustar_UB = repmat(u_UB,Nh,1);
zlim_LB = repmat(lift(FOM.C*FOM.x_LB,1),Nh+1,1);
zlim_UB = repmat(lift(FOM.C*FOM.x_UB,1),Nh+1,1);
for t = 1:tf/dt
    % Run finite-horizon optimization
    z0 = lift(y,t);
    cvx_begin quiet
    cvx_precision low
    variables z((Nh+1)*N,1) ustar(Nh*m,1)
    minimize (quad_form(z,G) + g'*z + quad_form(ustar,H) + h'*ustar)
    subject to
        % Dynamics
        z(N+1:end) == Ahat_blk*z(1:end-N) + Bhat_blk*ustar;
        
        % Initial conditions
        z(1:N) == z0;
        
        % Boundary conditions
        ustar_LB <= ustar; ustar <= ustar_UB;
        zlim_LB <= z; z <= zlim_UB;  % INFEASIBLE
    cvx_end
    
    % Fucker's at it again
    if isnan(ustar)
        ustar(1:m) = zeros(m,1);
    end
    
    % Update control sequence
    u = [u ustar(1:m)];
    
    % Progress dynamics
    [x_next, y_next] = full_dynamics(x(:,end),u(:,end));
    x = [x x_next];
    y = [y y_next];
end
end




