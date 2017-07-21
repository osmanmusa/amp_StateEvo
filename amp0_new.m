function [x, tau2, mse] = amp0_new(A, y, x_true, inter_max)

[M, N] = size(A);
delta = M/N;

tau2 = zeros(inter_max, 1);
mse = zeros(inter_max, 1);

z(:,1) = y;
x(:,1) = zeros(N,1);
tau2(1) = 1/M*(norm(z(:,1),2)^2);
x(:,2) = wthresh(A'*z(:,1) + x(:,1),'s',sqrt(tau2(1)));

mse(1) = norm(x_true - x,2).^2/N;

for t=2:inter_max+1

    eta_prime = 1/N * nnz(wthresh(A'*z(:,t-1) + x(:,t),'s',sqrt(tau2(t-1))));
    tau2(t) = tau2(t-1)/delta * eta_prime;
    
    
    eta_prime = 1/N * nnz(wthresh(A'*z(:,t-1) + x(:,t-1),'s',sqrt(tau2(t-1))));
    z(:,t) = y - A*x(:,t) + 1/delta * z(:,t-1) * eta_prime;

%     tau2(t) = 1/M*(norm(z(:,t),2)^2);
    
    x(:,t+1) = wthresh(A'*z(:,t) + x(:,t),'s',sqrt(tau2(t)));
    
    mse(t) = norm(x_true - x(:,t),2).^2/N;
      
    fprintf('AMP MSE = %f \n', 10*log10(mse(t)));   
end

fprintf('AMP final MSE = %f ', 10*log10(norm(x_true - x(:,inter_max+2),2).^2/N));
