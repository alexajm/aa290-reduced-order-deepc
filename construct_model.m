% Desired dimensions
n = 30;
m = 1;
p = 10;

% Limitations
max = 1;
min = -max;

% Randomize process dynamics
A = zeros(n,n);
while rank(A) ~= n
    A = (max - min)*rand(n,n) + min;
end

% a1 = 2;
% a2 = 1;
% A = diag(a1*ones(1,n)) + diag(a2*ones(1,n-1),1) + diag(a2*ones(1,n-1),-1);

% Pre-process matrix powers
Apow = cell(n,1);
for i = 0:n-1
    Apow{i+1} = A^i;
end

% Randomize control dynamics
% B = zeros(n,m);
% while ~is_controllable(Apow,B)
%     B = (max - min)*rand(n,m) + min;
% end

B = eye(n,m);

% Randomize observation dynamics
C = zeros(p,n);
while ~is_observable(Apow,C)
    C = (max - min)*rand(p,n) + min;
end

% Control limitations
u_UB = ones(p,1);
u_LB = -u_UB;

% Save model
FOM.A = A;
FOM.B = B;
FOM.C = C;
name = num2str(n) + "-" + num2str(m) + "-" + num2str(p) + "-synthetic.mat";
save(name,'FOM','u_UB','u_LB');

% Controllability
function cont = is_controllable(Apow,B)
    n = size(B,1);
    Cn = [];
    for i = 1:n
        Cn = [Cn, Apow{i} * B];
    end
    cont = rank(Cn) == n;
end

% Observability
function obsv = is_observable(Apow,C)
    n = size(C,2);
    Ob = [];
    for i = 1:n
        Ob = [Ob; C * Apow{i}];
    end
    obsv = rank(Ob) == n;
end




