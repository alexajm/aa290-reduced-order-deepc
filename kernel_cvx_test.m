m = 1;
p = 1;
lag = 2;
w = [1 2 3 4 5 6 7 8];

cvx_begin
variables R(1,(m+p)*lag)
minimize sum(R(end-(m+p):end))
subject to
    for j = 0:2
        0 == R * w((m+p)*j+1:(m+p)*(j+lag));
    end
cvx_end