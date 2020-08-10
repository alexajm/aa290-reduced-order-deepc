%% Setup

% controllers = ["deepc", "robust_deepc", "mpc", "koopman", "pgadp", "lqr"];
% files = ["large_synthetic", "distillation_column", "36-5-10-synthetic"];

global process_cov measure_cov FOM
file = 'large_synthetic';
tf = 8;                    % time to simulate until
dt = .2;                    % length of timestep
tspan = 0:dt:tf;
process_cov = 0;            % covariance of process noise
measure_cov = 0;            % covariance of measurement noise

% Load model
load(strcat(file, '.mat'));
n = size(FOM.A,1);       % state dimensionality
m = size(FOM.B,2);       % number of inputs
p = size(FOM.C,1);       % number of outputs

% Trajectory, initial conditions, noise settings
traj = zeros(n*(tf/dt + 1),1);
x0 = .5 * ones(n,1);
% x0 = (FOM.x_UB - FOM.x_LB).*rand(n,1) + FOM.x_LB;
process_cov = process_cov * eye(n);
measure_cov = measure_cov * eye(p);


%% DeePC vs MPC vs LQR
% Compare DeePC, MPC, and LQR control of a given system.

% Run DeePC
Nf = 10;        % forward time window (future)
Np = n;         % backward time window (past)
[H_full, ~, ~, ~] = deepc_offline(file,Np,Nf);
figure(); hold on
for H = [H_full]
    % Run
    [x,~] = deepc(x0, traj, tf, dt, H, Np, Nf);

    % Output results
    dist = zeros(1,size(x,2));
    for i = 1:length(dist)
        dist(i) = norm(x(:,i) - traj(n*(i-1)+1:n*i));
    end
    plot(tspan,dist);
end

% Run MPC
% [xm,~] = mpc(x0, traj, tf, dt);
% dist = zeros(1,size(xm,2));
% for i = 1:length(dist)
%     dist(i) = norm(xm(:,i) - traj(n*(i-1)+1:n*i));
% end
% plot(tspan,dist);

% Run LQR
[xl,~] = lqr(x(:,1), traj, tf, dt);
dist = zeros(1,size(xl,2));
for i = 1:length(dist)
    dist(i) = norm(xl(:,i) - traj(n*(i-1)+1:n*i));
end
plot(tspan,dist);

% Plot settings
title(['regulation of ',file]);
xlabel('time');
ylabel('norm(x - r)');
ylim([0 max(dist)+1]);
legend('deepc (full)','lqr');

%% Aggregate errors over many DeePC runs with many sample sizes
% Run DeePC for many datasets and view how cost changes with respect to
% how many dimensions the state space is assumed to occupy.

% Settings
est_dimensions = n-5:n;  % estimated dimensionality of the state space
num_dims = length(est_dimensions);
num_runs = 1;           % number of datasets to evaluate over
Nf = 10;                 % forward time window (future)
Np = n;                  % backward time window (past)
Q = 1e3 * ones(p,1);     % cost matrix

% Calculate costs
costs = zeros(num_runs,num_dims);
for i = 1:num_runs
    [~, ~, ~, datasets] = deepc_offline(file,Np,Nf,est_dimensions);
    for j = 1:num_dims
        [x,~] = deepc(x0, traj, tf, dt, datasets{j}, Np, Nf);
        for k = 1:size(x,2)
            y = FOM.C * x(:,k);
            costs(i,j) = costs(i,j) + y'*Q*y;
        end
    end
end
avgs = mean(costs,1);
devs = std(costs,0,1);

% Calculate baseline cost (MPC or LQR)
[x,~] = lqr(x0, traj, tf, dt);
base = 0;
for k = 1:size(x,2)
    y = FOM.C * x(:,k);
    base = base + y'*Q*y;
end

% Plot
figure(); hold on
%errorbar(est_dimensions,avgs,devs);
plot([0, max(est_dimensions)],[base, base], '--')
title('deepc cost vs estimated dimension')
xlabel('estimated dimension')
ylabel('cost')
%xlim([min(est_dimensions) - .2, max(est_dimensions) + .2])

%% Simulate known control policy with Hankel dynamics
% Given an LQR control policy, simulate the system dynamics exclusively
% using the Hankel of the signal and compare to actual dynamics.

% Collect data
[H_full, ~, H_redux, ~] = deepc_offline(file,Np,Nf);
Up = H_full.Up;
Uf = H_full.Uf;
Yp = H_full.Yp;
Yf = H_full.Yf;

% Collect Np measurements prior to controller activation
x_minus = 2 * ones(n,1);
u_minus = zeros(m,1);
y_minus = [];
for t = -Np:0
    [x_minus, y_next] = full_dynamics(x_minus,u_minus);
    y_minus = [y_minus y_next];
end
x0 = x_minus(:,end);

% Run LQR
[x,u] = lqr(x0, traj, tf, dt);

% Construct actual dynamics
y = y_minus;
for t = 2:tf/dt
    y_next = FOM.C * x(:,t);
    y = [y y_next];
end

% Construct full order Hankel dynamics
yh = y_minus;
u = [zeros(p,Np) u];
for t = 1:(tf/dt - Nf)
    % Collect past data, control policy
    t0 = t+Np;
    u0 = reshape(u(:,t:t0-1),numel(u(:,t:t0-1)),1);
    y0 = reshape(yh(:,end-Np+1:end),numel(yh(:,end-Np+1:end)),1);
    uf = reshape(u(:,t0:t0+Nf-1),numel(u(:,t0:t0+Nf-1)),1);
    
    % Progress dynamics
    yh_next = hankel_dynamics(Up,Uf,u0,uf,Yp,Yf,y0);
    if t ~= tf/dt - Nf
        yh_next = yh_next(1:p);
    else
        yh_next = reshape(yh_next,p,Nf);
    end
    yh = [yh yh_next];
end

% Construct reduced order Hankel dynamics
Up = H_redux.Up;
Uf = H_redux.Uf;
Yp = H_redux.Yp;
Yf = H_redux.Yf;
yh_redux = y_minus;
for t = 1:(tf/dt - Nf)
    % Collect past data, control policy
    t0 = t+Np;
    u0 = reshape(u(:,t:t0-1),numel(u(:,t:t0-1)),1);
    y0 = reshape(yh_redux(:,end-Np+1:end),numel(yh_redux(:,end-Np+1:end)),1);
    uf = reshape(u(:,t0:t0+Nf-1),numel(u(:,t0:t0+Nf-1)),1);
    
    % Progress dynamics
    yh_redux_next = hankel_dynamics(Up,Uf,u0,uf,Yp,Yf,y0);
    if t ~= tf/dt - Nf
        yh_redux_next = yh_redux_next(1:p);
    else
        yh_redux_next = reshape(yh_redux_next,p,Nf);
    end
    yh_redux = [yh_redux yh_redux_next];
end

% Plot
figure(); hold on
tspan = (1-Np)*dt:dt:tf-dt;
plot(tspan,y(1:end-1))
plot(tspan,yh(1:end-1),'--')
plot(tspan,yh_redux(1:end-1),'-.')
title('control simulation for state-space and Hankel dynamics')
xlabel('time')
ylabel('output')
legend('actual','hankel','redux')
xlim([0 tf]);




