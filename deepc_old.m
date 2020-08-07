function [x, u] = deepc_old(x0, traj, tf, dt, H)
% Globals
global ROM FOM u_UB u_LB
global measure_cov

% Setup
Nf = 10;            % forward time window
Np = 1;             % backward time window
L = Np + Nf;
m = size(FOM.B,2);  % number of inputs
n = size(FOM.A,2);  % state dimensionality
p = size(FOM.C,1);  % number of outputs
lam_g = 1e2;
lam_y = 1e5;

% Collect input/output Hankels
Up = H.Up;
Uf = H.Uf;
Yp = H.Yp;
Yf = H.Yf;

% Construct performance matrix
Qbase = 1e2 * ones(p,1);
Qdiag = repmat(Qbase,Nf,1);
Q = diag(Qdiag);

% Track trajectory with DeePC
u = [];
x = x0;
u_in = [];
y_out = FOM.C*x;
ustar_LB = repmat(u_LB,Nf,1);
ustar_UB = repmat(u_UB,Nf,1);
for t = 1:tf/dt
    % Get reference trajectory for this window
    ref = traj(p*(t-1)+1:p*(t+Nf-1));
    
    % CVX w/ Kalman
    ustar = 0;
    if t > Np && kalman
        cvx_begin quiet
        cvx_precision default
        variables g(size(Up,2)+q,1) ustar(Nf*m,1)
        variables y(Nf*p,1) sig(length(y_out),1)
        minimize (lam_g*norm(g,1) + lam_y*norm(sig,1) +...
                  quad_form(y - ref, Q) +...
                  quad_form(ustar, eye(length(ustar))))
        subject to
            [[Up(1:length(u_in),:); Yp(1:length(y_out),:); Uf; Yf],mu] * g ==...
                [u_in; y_out; ustar; y] +...
                [zeros(size(u_in)); sig; zeros(Nf*(m+p),1)];
            ustar_LB <= ustar; ustar <= ustar_UB;
            %ylim_LB <= y; y <= ylim_UB;
            sig >= 0;
        cvx_end
        ustar = [Uf,mu(1:Nf*m)] * g;
    else
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
            ustar_LB <= ustar; ustar <= ustar_UB;
            %ylim_LB <= y; y <= ylim_UB;
            sig >= 0;
        cvx_end
        ustar = Uf * g;
    end
    
    % CVX is back on its shit
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
    
    % Kalman update
    if t > Np && kalman
        mu1 = mu(end-Nf+1:end-Nf+p);
        C = [Yf(1:p,:),mu1];
        K = (S*C') * inv(C*S*C' + measure_cov*eye(size(C,1)));
        M = K*(y_next - C*g);
        y_next - C*g;
        mu1 = mu1 + M(end-p+1:end);
        S = (eye(size(S)) - K*C)*S;
        mu(end-Nf+1:end-Nf+p) = mu1;
    end
end
end

% Nearest neighbor estimate of Shannon entropy
function E = entropy(H)
    n = size(H,2);
    E = 0;
    for i = 1:n
        % Get nearest neighbor distance rho
        wi = H(:,i);
        wj = H(:,1);
        rho = Inf;
        for j = 1:size(H,2)
            if i == j; continue; end
            if norm(wi - H(:,j),2) < rho
                wj = H(:,j);
                rho = norm(wi - wj,2);
            end
        end
        
        % Add to entropy
        E = E + log(n*rho);
    end
    E = E / n;
end

function avg = avg_delta(deltas,w)
    avg = zeros(length(deltas)-w,1);
    for i = 1:w
        avg = avg + deltas(i:end-w+i-1);
    end
    avg = avg ./ w;
end

function lam = lambda(b)
    c = 2*(b + 1);
    d = 8*b;
    e = (b+1) + sqrt(b^2 + 14*b + 1);
    lam = sqrt(c + d / e);
end




