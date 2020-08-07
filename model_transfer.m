% load('tubular_reactor_FOM.mat')
% load('tubular_reactor_ROM.mat')
% load('large_synthetic_FOM.mat')
% load('large_synthetic_ROM.mat')
load('distillation_column_FOM.mat')
load('distillation_column_ROM.mat')
FOM.A = FOM.Af;
FOM.B = FOM.Bf;
FOM.Bw = FOM.Bfw;
FOM.C = FOM.Cf;
FOM.H = FOM.Hf;
FOM = rmfield(FOM,'Af');
FOM = rmfield(FOM,'Bf');
FOM = rmfield(FOM,'Bfw');
FOM = rmfield(FOM,'Cf');
FOM = rmfield(FOM,'Hf');

% large synthetic
% o = size(FOM.H,1);
% m = size(FOM.B, 2);
% ROM.x_UB = ones(1,o); ROM.x_UB = ROM.x_UB';
% FOM.x_UB = ones(1,size(FOM.A,1)); FOM.x_UB = FOM.x_UB';
% ROM.x_LB = -ROM.x_UB;
% FOM.x_LB = -FOM.x_UB;
% u_UB = 0.25*ones(m,1);
% u_LB = -u_UB;

% tubular reactor
% P.N = 300; % spatial discretization parameter
% P.v = 0.1;
% P.L = 1;
% P.dz = P.L/P.N;
% P.Tnorm = 340;
% P.Cnorm = 0.02;
% P.k0 = 10^6;
% P.E = 11250;
% P.R = 1.986;
% P.Cin = 0.02;
% P.Tin = 340;
% P.Gr = 4.25*10^9;
% P.Hr = 0.2;
% TJ1 = 374.6;
% TJ2 = 310.1;
% TJ3 = 325.2;
% [Cstar, Tstar] = reactorSteadyState(TJ1, TJ2, TJ3, P);
% Tmax = 395;
% Tmin = 300;
% Tnorm = 340;
% umax = 395;
% umin = 300;
% FOM.x_UB = (Tmax - Tstar)/Tnorm; FOM.x_UB = FOM.x_UB';
% FOM.x_LB = (Tmin - Tstar)/Tnorm; FOM.x_LB = FOM.x_LB';
% ROM.x_UB = FOM.x_UB(1:30:end);
% ROM.x_LB = FOM.x_LB(1:30:end);
% u_UB = [(umax - TJ1)/Tnorm, (umax - TJ2)/Tnorm, (umax - TJ3)/Tnorm]';
% u_LB = [(umin - TJ1)/Tnorm, (umin - TJ2)/Tnorm, (umin - TJ3)/Tnorm]';
% FOM.x_UB = FOM.H * FOM.x_UB;
% FOM.x_LB = FOM.H * FOM.x_LB;
% ROM.x_UB = ROM.H * ROM.x_UB;
% ROM.x_LB = ROM.H * ROM.x_LB;

% distillation column
zc = [0.01, 0.01, 1, 1]; % constraints on real z = [y_D, x_B, M_D, M_B]
uc = [1, 1, 1, 1]; % constraints on real u = [L, V, D, B]
duc = [0.2, 0.2, 0.2, 0.2]; % constraints on change in u
FOM.x_UB = [zc, uc]';
FOM.x_LB = -FOM.x_UB;
u_UB = duc';
u_LB = -u_UB;

%save('tubular_reactor.mat','ROM','FOM','u_UB','u_LB');
%save('large_synthetic.mat','ROM','FOM','u_UB','u_LB');
save('distillation_column.mat','ROM','FOM','u_UB','u_LB');




