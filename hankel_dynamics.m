function y = hankel_dynamics(Up,Uf,u0,u,Yp,Yf,y0)
    A = [Up(1:length(u0),:); Uf; Yp(1:length(y0),:)];
    b = [u0; u; y0];
    g = pinv(A) * b;
    y = Yf * g;
end