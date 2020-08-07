function H = hankel_dp(H0,u,y)
    % Retrieve parameters
    m = size(u,1);
    p = size(y,1);
    L = size(H0,1) / (m + p);
    T = size(H0,2) + L - 1;
    
    
end