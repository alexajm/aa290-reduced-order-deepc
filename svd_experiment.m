L = 2;
u = [1 2 3]';
H3 = hankel(L,3,u);
[~,S3,~] = svd(H3);

u = [1 2 3 4]';
H4 = hankel(L,4,u);
[~,S4,~] = svd(H4);

u = [1 2 3 4 5]';
H5 = hankel(L,5,u);
[~,S5,~] = svd(H5);

norm(S3,'fro')
norm(S4,'fro')
norm(S5,'fro')

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