function x = robust_deepc(x0, traj, tf, dt)
% Globals
global ROM FOM u_UB u_LB

% Setup
% Nf is length of the forward-looking time window
% Np is length of the backward-looking time window
Nf = 10;
Np = 1;
m = size(FOM.B,2);
n = size(FOM.A,2);
p = size(FOM.C,1);
lam_g = 1;
lam_y = 1e5;
L = Np + Nf;
ulim = 100;
ylim = 2;
ep = 1e-3;
N = 10;
alpha = .1;
r = 2;
q = 1 / (1 - (1/r));

% Get persistently excited T
% T represents amount of offline data the model is fit to
T = L * ((m*L + 1)*(n + 1) - 1);
%T = 50;

% Collect offline data
u_off = 2 .* ulim .* rand(m,T) - ulim;
y_off = cell(N,1);
for i = 1:N
    x = zeros(n,1);
    %yi = FOM.C * x;
    yi = [];
    for j = 1:T
        % Apply dynamics
        [x_next, yi_next] = full_dynamics(x, u_off(:,j));
        
        % Save results
        yi = [yi; yi_next];
        x = x_next;
    end
    y_off{i} = yi;
end

% Process offline data
u_off = reshape(u_off,size(u_off,1)*size(u_off,2),1);
Pu = page(L,T,u_off);
Up = Pu(1:Np*m,:);
Uf = Pu(Np*m+1:end,:);
Yp = cell(N,1);
Yf = cell(N,1);
for i = 1:N
    Py = page(L,T,y_off{i});
    Yp{i} = Py(1:Np*p,:);
    Yf{i} = Py(Np*p+1:end,:);
end

% Cost matrices
Q = eye(size(Uf,1));

% Get Lipschitz constants
% NOTE: these are guesses, not analytical calculations
Lobj = 1e2;
Lcon = 1e2;

% Track trajectory
u = [];
x = x0;
u_in = u_off(end-Np*m+1:end);
y_out = [y_off{N}(end-(Np-1)*p+1:end); FOM.C*x];
% ylim_LB = repmat(FOM.C*FOM.x_LB,Nf,1);
% ylim_UB = repmat(FOM.C*FOM.x_UB,Nf,1);
ylim_LB = repmat(FOM.C*pinv(FOM.H)*FOM.x_LB,Nf,1);
ylim_UB = repmat(FOM.C*pinv(FOM.H)*FOM.x_UB,Nf,1);
for t = 1:tf/dt
    % Get reference trajectory for this window
    ref = traj(p*(t-1)+1:p*(t+Nf-1));
    
    % Run CVX
    % NOTE: f1, f2 don't follow paper precisely
    cvx_begin
    cvx_precision low
    variables g(size(Up,2),1) tau s(N,1) sig(size(u_in))
    expression f1
        f1 = quad_form(Uf*g,Q);
    expression f2(N)
        for i = 1:N
            f2(i) = 1e3 * norm(Yf{i}*g - ref, 2);
        end
    expression f3(N)
        for i = 1:N
            f3(i) = 1e3 * norm(Yp{i}*g - y_out, 2);
        end
    expression h(N)
        for i = 1:N
            h(i) = max([Yf{i}*g; -Yf{i}*g] - [ylim_UB; -ylim_LB]);
        end
    minimize (f1 + sum(f2)/N + sum(f3)/N +...
        Lobj*ep*norm(g,q) + 1e10*norm(sig,2));
    subject to
        Up * g == u_in + sig;
        -ulim <= Uf*g; Uf*g <= ulim;
        0 >= -tau*alpha + Lcon*ep*norm(g,q) + norm(s,1)/N;
        for i = 1:N
            tau + h(i) <= s(i);
        end
        s >= 0;
    cvx_end

    % Update control sequence
    ustar = Uf * g;
    if isnan(ustar)
        ustar(1:m) = zeros(m,1);
    end
    u = [u ustar(1:m)];
    
    % Progress dynamics
    [x_next, y_next] = full_dynamics(x(:,end),u(:,end));
    x = [x x_next];

    % Update input/output data
    u_in = [u_in(Np*m+1:end); u(:,end)];
    y_out = [y_out(Np*p+1:end); y_next];
    if isnan(y_out)
        y_out
    end
end
end

function P = page(L,T,d)
    n = length(d) / T;
    rows = L * n;
    cols = floor(T/L);
    P = zeros(rows,cols);
    for i = 0:L-1
        for j = 0:cols-1
            p1 = n*i + 1;
            p2 = n*(i + 1);
            d1 = n*(i + L*j) + 1;
            d2 = n*(i + L*j + 1);
            P(p1:p2,j+1) = d(d1:d2);
        end
    end
end

function ksi = get_ksi(Yp,Yf)
    Ymat = [Yp; Yf];
    ksi = [];
    for i = 1:size(Ymat,1)
        ksi = [ksi Ymat(i,:)'];
    end
    ksi = ksi';
end




