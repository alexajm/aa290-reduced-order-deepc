function [x,u] = mpc(x0, traj, tf, dt)
% Globals
global FOM u_UB u_LB

% Settings
N = 10;
A = FOM.A;
B = FOM.B;
C = FOM.C;
m = size(B,2);
n = size(A,2);
p = size(C,1);

% Block-ify dynamics
A_blk = [];
B_blk = [];
C_blk = [];
for i = 1:N
    A_blk = blkdiag(A_blk,A);
    B_blk = blkdiag(B_blk,B);
    C_blk = blkdiag(C_blk,C);
end

% Construct cost matrices
Qbase = 100 * ones(p,1);
Qdiag = repmat(Qbase,N,1);
Q = diag(Qdiag);
Rbase = 10 * ones(m,1);
Rdiag = repmat(Rbase,N,1);
R = diag(Rdiag);

u = [];
x = x0;
mu = x0;
S = .1*eye(n);
ustar_LB = repmat(u_LB,N,1);
ustar_UB = repmat(u_UB,N,1);
% ylim_LB = repmat(FOM.C*FOM.x_LB,N,1);
% ylim_UB = repmat(FOM.C*FOM.x_UB,N,1);
for t = 1:tf/dt
    % Shorten horizon as you approach end of trajectory
    while n*(t+N) > length(traj)
        N = N - 1;
        A_blk = A_blk(1:N*n,1:N*n);
        B_blk = B_blk(1:N*n,1:N*m);
        C_blk = C_blk(1:N*p,1:N*n);
        Q = Q(1:N*p,1:N*p);
        R = R(1:N*m,1:N*m);
        ustar_LB = repmat(u_LB,N,1);
        ustar_UB = repmat(u_UB,N,1);
    end
    
    % Get reference trajectory for this window
    ref = C_blk * traj(n*t+1:n*(t+N));
    
    % Run CVX
    cvx_begin quiet
    cvx_precision default
    variables ustar(N*m,1) z((N+1)*n,1) y(N*p,1)
    minimize (quad_form(y - ref,Q) + quad_form(ustar,R))
    subject to
        % Dynamics
        z(n+1:end) == A_blk*z(1:end-n) + B_blk*ustar;
        y == C_blk * z(1:end-n);
        
        % Initial conditions
        z(1:n) == mu;
        
        % Boundary conditions
        %ustar_LB <= ustar; ustar <= ustar_UB;
        %ylim_LB <= y; y <= ylim_UB;
    cvx_end
    
    % Fix this shite
    % See if this works, tf can't CVX do it ://
    if isnan(ustar)
        ustar(1:m) = zeros(m,1);
    end
    
    % Update control sequence
    u = [u ustar(1:m)];
    
    % Progress dynamics
    [x_next, y_next] = full_dynamics(x(:,end),u(:,end));
    x = [x x_next];
    
    % Predict, update state estimate
    % NOTE: full state knowledge assumed
%     [mu_prior, S_prior] = predict(mu, S, u(:,end), A);
%     [mu, S] = update(mu_prior, S_prior, C, y_next);
    mu = x_next;
    
    if isnan(mu)
        mu
    end
end
end

% Kalman prediction step
function [mu_prior, S_prior] = predict(mu_prev,S_prev,u,A)
    global process_cov
    mu_prior = full_dynamics(mu_prev,u);
    S_prior = A*S_prev*A' + process_cov*eye(size(S_prev));
end

% Kalman update step
function [mu_post, S_post] = update(mu_prior,S_prior,C,y)
    global measure_cov
    K = (S_prior*C') * inv(C*S_prior*C' + measure_cov*eye(size(C,1)));
    mu_post = mu_prior + K*(y - C*mu_prior);
    S_post = (eye(size(S_prior)) - K*C)*S_prior;
end




