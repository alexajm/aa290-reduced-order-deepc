%% DeePC vs MPC vs LQR
% Compare DeePC, MPC, and LQR control of a given system.

% controllers = ["deepc", "robust_deepc", "mpc", "koopman", "pgadp", "lqr"];
% files = ["large_synthetic", "distillation_column", "36-5-10-synthetic"];

% Control, model settings
global process_cov measure_cov FOM
file = 'large_synthetic';
tf = 10;                    % time to simulate until
dt = .1;                    % length of timestep
tspan = 0:dt:tf;
process_cov = 0;            % covariance of process noise
measure_cov = 0;            % covariance of measurement noise

% Load model
load(strcat(file, '.mat'));
n = size(FOM.A,1);
p = size(FOM.C,1);

% Trajectory, initial conditions, noise settings
traj = zeros(n*(tf/dt + 1),1);
x0 = .5 * ones(n,1);
% x0 = (FOM.x_UB - FOM.x_LB).*rand(n,1) + FOM.x_LB;
process_cov = process_cov * eye(n);
measure_cov = measure_cov * eye(p);

% Run DeePC
Nf = 10;        % forward time window (future)
Np = 1;         % backward time window (past)
[H_full, H_redux, H_det] = deepc_offline(file,Np,Nf);
figure(); hold on
for H = [H_full, H_redux, H_det]
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
[xm,~] = mpc(x0, traj, tf, dt);
dist = zeros(1,size(xm,2));
for i = 1:length(dist)
    dist(i) = norm(xm(:,i) - traj(n*(i-1)+1:n*i));
end
plot(tspan,dist);

% Run LQR
[xl,~] = lqr(x0, traj, tf, dt);
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
legend('deepc (full)','deepc (redux)','deepc (det)','mpc','lqr');

%% Aggregate errors over many DeePC runs with many sample sizes
% Run DeePC for many datasets and view how error changes with respect to
% the proportion of data used

% Settings
global FOM
sample_sizes = .4:.1:1;  % proportion of persistently exciting data to use
num_sizes = length(sample_sizes);
num_runs = 10;           % number of datasets to evaluate over
Nf = 10;                 % forward time window (future)
Np = 1;                  % backward time window (past)
file = 'large_synthetic';
load(strcat(file, '.mat'));
n = size(FOM.A,1);       % state dimensionality
m = size(FOM.B,2);       % number of inputs
p = size(FOM.C,1);       % number of outputs
x0 = .5 * ones(n,1);     % initial state
traj = zeros(n*(tf/dt + 1),1); % reference trajectory
tf = 8;                  % time to simulate until
dt = .2;                 % length of timestep

% Calculate errors
errors = zeros(num_runs,num_sizes);
for i = 1:num_runs
    [~, ~, ~, datasets] = deepc_offline(file,Np,Nf,sample_sizes);
    for j = 1:num_sizes
        x = deepc(x0, traj, tf, dt, datasets{j}, Np, Nf);
        for k = 1:size(x,2)
            errors(i,j) = errors(i,j) + norm(x(:,k));
        end
    end
end
avgs = mean(errors,1);
devs = std(errors,0,1);

% Plot
figure()
errorbar(sample_sizes,avgs,devs);
title('deepc error vs sample size')
xlabel('proportion of persistently exciting data used')
ylabel('error')
xlim([min(sample_sizes) - .2, max(sample_sizes) + .2])

%% Compare real dynamics, Hankel dynamics
% NOTE: not currently compatible with recent updates to code

% Construct full order Hankel dynamics
Np = 1;
Nf = 10;
[H_full, H_redux] = deepc_offline(file,Np,Nf);
Up = H_full.Up;
Uf = H_full.Uf;
Yp = H_full.Yp;
Yf = H_full.Yf;
yh = FOM.C * x0;
y = FOM.C * x0; 
for t = 2:(tf/dt - Nf)
    % Update Hankel measurement
    up = [];
    if t-1 < Np
        up = reshape(u(:,1:t),numel(u(:,1:t)),1);
    else
        up = reshape(u(:,t-Np:t-1),numel(u(:,t-Np:t-1)),1);
    end
    uf = reshape(u(:,t:t+Nf-1),numel(u(:,t:t+Nf-1)),1);
    y0 = [];
    if t <= Np
        y0 = reshape(y,numel(y),1);
    else
        y0 = reshape(y(:,end-Np+1:end),numel(y(:,end-Np+1:end)),1);
    end
    yh_next = hankel_dynamics(Up,Uf,up,uf,Yp,Yf,y0);
    yh_next = yh_next(1:p);
    yh = [yh yh_next];

    % Update actual measurement
    y_next = FOM.C * x(:,t);
    y = [y y_next];
end

% Construct reduced order Hankel dynamics
Up = H_full.Up;
Uf = H_full.Uf;
Yp = H_full.Yp;
Yf = H_full.Yf;
yh_low = FOM.C * x0;
for t = 2:(tf/dt - Nf)
    % Update Hankel measurement
    up = [];
    if t-1 < Np
        up = reshape(u(:,1:t),numel(u(:,1:t)),1);
    else
        up = reshape(u(:,t-Np:t-1),numel(u(:,t-Np:t-1)),1);
    end
    uf = reshape(u(:,t:t+Nf-1),numel(u(:,t:t+Nf-1)),1);
    y0 = [];
    if t <= Np
        y0 = reshape(y,numel(y),1);
    else
        y0 = reshape(y(:,end-Np+1:end),numel(y(:,end-Np+1:end)),1);
    end
    yh_low_next = hankel_dynamics(Up,Uf,up,uf,Yp,Yf,y0);
    yh_low_next = yh_low_next(1:p);
    yh_low = [yh_low yh_low_next];
end
figure()
dist = zeros(1,size(y,2));
dist_h = zeros(1,size(yh,2));
dist_low = zeros(1,size(yh_low,2));
for i = 1:length(dist)
    ob = FOM.C*traj(n*(i-1)+1:n*i);
    dist(i) = norm(y(:,i) - ob);
    dist_h(i) = norm(yh(:,i) - ob);
    dist_low(i) = norm(yh_low(:,i) - ob);
end
full_error = norm(dist_h - dist);
redux_error = norm(dist_low - dist);
disp(['full error ',num2str(full_error),', reduced error ',...
      num2str(redux_error),' over ',num2str(tf/dt),' timesteps']);

% Plot
T = (0:length(dist)-1)*dt;
plot(T,dist,T,dist_h,T,dist_low)
legend('actual','full hankel','redux hankel');
title(['application of ',control_str,' policy to Hankel for ',file]);
xlabel('time');
ylabel('norm(y - r)');
%ylim([0 max([dist; dist_h])+1]);







