function experiment(file, control_type, tf, dt, x0)
% Setup
% file = 'large_synthetic';
% control_type = 'robust_deepc';
% tf = 20;
% dt = 1;

% Load full and reduced-order models
global FOM ROM
load(strcat(file, '_FOM.mat'));
load(strcat(file, '_ROM.mat'));

% Select desired trajectory (for now: regulation)
traj = zeros(size(ROM.A,1)*tf,1);

% Collect controller
controller = collect_controller(control_type);

% Run simulation
x = controller(x0, traj, tf, dt);

% Output results
figure()
T = 0:dt:tf;
dist = zeros(1,size(x,2));
for i = 1:length(dist)
    dist(i) = norm(x(:,i));
end
plot(T,dist)
title(['regulation of ',file,' over time using ',control_type,' controller']);
xlabel('time');
ylabel('norm(x - r)');
ylim([0 max(dist)+1]);

% Retrieve handle for controller function
function controller = collect_controller(control_type)
switch control_type
    case 'deepc'
        controller = @deepc;
    case 'koopman'
        controller = @koopman;
    case 'pgadp'
        controller = @pgadp;
    case 'robust_deepc'
        controller = @robust_deepc;
    otherwise
        error('Invalid controller (control_type = %s)',control_type);
end
end
end




