% DeePC offline data collection
%
% inputs
%       -file: file to access model from
%       -Np: backward time window (past)
%       -Nf: forward time window (future)
%       -sample_sizes: proportions of full data to generate Hankels for
%
% outputs (structs of input/output Hankels Up, Yp, Uf, Yf)
%       -H_full: persistently exciting data
%       -H_redux: data reduced by parallel analysis
%       -H_det: data necessary to precisely determine DeePC optimization
%       -datasets: cell of Hankels corresponding to proportions of full
%           data described by sample_sizes
%
function [H_full,H_redux,H_det,datasets] = deepc_offline(file,Np,Nf,sample_sizes)
% Globals
global FOM

% Setup
if nargin == 3
    sample_sizes = [];
end

% Important values
load(strcat(file, '.mat'));
m = size(FOM.B,2);
n = size(FOM.A,2);
p = size(FOM.C,1);
L = Np + Nf;

% Input constraints
% NOTE: currently hardcoded to large_synthetic
u_UB = .25;
u_LB = -.25;

% Get persistently excited sample size T
TFOM = (m + 1)*(Np + Nf + n) - 1;
T_full = TFOM;
T_redux = -1;
T_det = Np*(m+1) + Nf*(m+p+1) - 1;

% Collect offline data
u_off = zeros(T_full*m,1);
y_off = zeros(T_full*p,1);
x = zeros(n,1);
w = L;
past_ranks = -ones(w,1);
%figure();
for t = 1:T_full
    % Apply dynamics
    u_next = (u_UB - u_LB).*rand(m,1) + u_LB;
    [x_next, y_next] = full_dynamics(x, u_next);

    % Save results
    u_off(m*(t-1)+1:m*t) = u_next;
    y_off(p*(t-1)+1:p*t) = y_next;
    x = x_next;

    % Get approximate Hankel rank
    if t > L && T_redux < 0
        % Construct Hankel of signal
        H = hankel(L,t,[u_off(1:t*m); y_off(1:t*p)]);

        % Parallel analysis
        [lat, lat_hi, ~] = pa_test(H',200);
        low_rank = sum(lat >= lat_hi');
        past_ranks = [past_ranks(2:end); low_rank];
        if all(past_ranks == past_ranks(end))
            disp(['window behavior is approximately ', num2str(low_rank), '-dimensional']);
            disp(['determined with ', num2str(t), ' samples (reduced from ', num2str(TFOM), ')']);
            T_redux = t;
            %break
        end
    end
end

% Process full offline data
u_off = u_off(1:T_full*m);
Hu = hankel(L,T_full,u_off);
H_full.Up = Hu(1:Np*m,:);
H_full.Uf = Hu(Np*m+1:end,:);
y_off = y_off(1:T_full*p);
Hy = hankel(L,T_full,y_off);
H_full.Yp = Hy(1:Np*p,:);
H_full.Yf = Hy(Np*p+1:end,:);

% Process reduced offline data
if T_redux == -1; T_redux = T_full; end
ur_off = u_off(1:T_redux*m);
Hu = hankel(L,T_redux,ur_off);
H_redux.Up = Hu(1:Np*m,:);
H_redux.Uf = Hu(Np*m+1:end,:);
yr_off = y_off(1:T_redux*p);
Hy = hankel(L,T_redux,yr_off);
H_redux.Yp = Hy(1:Np*p,:);
H_redux.Yf = Hy(Np*p+1:end,:);

% Process deterministic offline data
ud_off = u_off(1:T_det*m);
Hu = hankel(L,T_det,ud_off);
H_det.Up = Hu(1:Np*m,:);
H_det.Uf = Hu(Np*m+1:end,:);
yd_off = y_off(1:T_det*p);
Hy = hankel(L,T_det,yd_off);
H_det.Yp = Hy(1:Np*p,:);
H_det.Yf = Hy(Np*p+1:end,:);

% Process full datasets
num_samples = length(sample_sizes);
datasets = cell(num_samples);
for i = 1:num_samples
    sze = floor(sample_sizes(i) * T_full);
    
    % Get u
    ui_off = u_off(1:sze*m);
    Hu = hankel(L,sze,ui_off);
    Hi.Up = Hu(1:Np*m,:);
    Hi.Uf = Hu(Np*m+1:end,:);
    
    % Get y
    yi_off = y_off(1:sze*p);
    Hy = hankel(L,sze,yi_off);
    Hi.Yp = Hy(1:Np*p,:);
    Hi.Yf = Hy(Np*p+1:end,:);
    datasets{i} = Hi;
end

disp(['sample size: ',num2str(T_full),' (full) ',num2str(T_redux),...
    ' (redux) ',num2str(T_det),' (det)'])




