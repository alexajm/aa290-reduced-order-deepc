% Construct Hankel matrix of order L from scratch given T samples of data
% contained in d
function H = hankel(L,T,d)
    n = length(d) / T;
    H = zeros(L*n,T - L + 1);
    for i = 0:L-1
        for j = 0:T - L
            h1 = n*i + 1;
            h2 = n*(i+1);
            d1 = n*(i+j) + 1;
            d2 = n*(i+j+1);
            H(h1:h2,j+1) = d(d1:d2);
        end
    end
end