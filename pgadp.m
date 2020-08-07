function x = pgadp(x0, traj, tf, dt)
% Globals
global ROM FOM u_LB u_UB

% Setup
alpha = .5;
num_samples = 1e3;
%num_samples = 5;
m = size(FOM.B,2);
n = size(FOM.A,2);
p = size(FOM.C,1);
T = tf / dt;

% Collect offline data
%S = zeros(num_samples, 3);
S = cell(num_samples,3);
x = zeros(n,1);
y = FOM.C * x;
for i = 1:num_samples
    % Generate sample
    u = (u_UB - u_LB).*rand(m,1) + u_LB;
    [x_next, y_next] = full_dynamics(x,u);
    S(i,:) = {y, u, y_next};
    
    % Update state
    x = x_next;
    y = y_next;
end

% Compute A0, y0, F0
A0 = 0;
z0 = zeros(size(weight_q(zeros(p,1),zeros(m,1))));
F0 = 0;
for i = 1:num_samples
    % Get sample
    xi = cell2mat(S(i,1));
    ui = cell2mat(S(i,2));
    
    % Update values
    A0 = A0 + weight_q(xi,ui)' * critic(xi,ui);
    z0 = z0 + cost(xi,ui) .* weight_q(xi,ui);
    F0 = F0 + weight_u(xi)' * actor(xi);
end

% Real-time actor-critic adaptive control
x = x0;
y = FOM.C * x;
gamma = zeros(size(actor(y),1),m);
%gamma = .1 * ones(size(actor(y),1),m);
u = zeros(m,1);
%u = .2;
s = cell(T,3);
for t = 1:T
    % Policy evaluation
    % Apply control
    [x_next, y_next] = full_dynamics(x(:,end),u);
    s(i,:) = {y, u, y_next};
    
    % Get new control policy
    u_next = (actor(y)' * gamma)';

    % Calculate important values
    A = A0 + weight_q(y,u)' * critic(y,u);
    B = weight_q(y,u)' * critic(y_next,u_next);
    for i = 1:num_samples
        yi = cell2mat(S(i,1));
        ui = cell2mat(S(i,2));
        yi_next = cell2mat(S(i,3));
        B = B + weight_q(yi,ui)' * critic(yi_next,(actor(yi_next)'*gamma)');
    end
    z = z0 + cost(y,u);
    %rcond(A)
    rcond(B);
    %rcond(A - B)
    beta = (A - B) \ z;
    
    % Policy improvement
    F = F0 + weight_u(y)' * actor(y);
    G = weight_u(y) * grad_critic(y,u_next)';
    for i = 1:num_samples
        yi = cell2mat(S(i,1));
        G = G + weight_u(yi)' * grad_critic(yi,(actor(yi)'*gamma)');
    end
    
%     size(weight_u(y))
%     size(actor(y))
%     size(grad_critic(y,actor(y)'*gamma))
%     size(F)
%     size(G)
    
    %rcond(F)
    %gamma
    gamma = gamma - alpha * (F \ G) * beta
    
    % Update state
    x = [x x_next];
    y = y_next;
    u = u_next;
end
end

% Positive-definite cost function
function J = cost(x, u)
global FOM u_LB u_UB
% S = (x' * x) + 1e0*(sum(x <= FOM.C*pinv(FOM.H)*FOM.x_LB) +...
%     sum(x >= FOM.C*pinv(FOM.H)*FOM.x_UB));
% W = (u' * u) + 1e10*(sum(u <= u_LB) + sum(u >= u_UB));
S = x' * x;
if any(x <= FOM.C*pinv(FOM.H)*FOM.x_LB) ||...
        any(x >= FOM.C*pinv(FOM.H)*FOM.x_UB)
    S = S * 1e3;
end
W = u' * u;
if any(u <= u_LB) || any(u >= u_UB)
    W = W * 1e3;
end
J = S + W;
end

% State-action weights
function Wq = weight_q(x, u)
Wq = [x; u; x*u; x.^2];
%Wq = [u.^2; u; x.^2; x];
end

% State weights
function Wu = weight_u(x)
global FOM
Wu = [x; x.^2; x.^3; x.^4];
%Wu = [x.^2; x.^4];
%Wu = repmat(Wu,1,size(FOM.B,2));
end

% Actor function vector
function psi = critic(x, u)
%psi = [x; x.^2; u; u.^2];
psi = [x; u; x*u; x.^2];
%psi = repmat(psi,1,size(u,1));
end

% Actor gradient vector
function grad_psi = grad_critic(x, u)
%grad_psi = [zeros(size(x)); x; zeros(size(x)); x.^2];
grad_psi = [zeros(size(x)); ones(size(u)); x; zeros(size(x))];
end

% Critic function vector
function phi = actor(x)
global FOM
phi = [x; x.^2; x.^3; x.^4];
%phi = [x; x.^3];
%phi = repmat(phi,1,size(FOM.B,2));
end




