function [x, ampsim_tau2, ampsim_MSE] = amp0(A, y, x_true, inter_max)

[M, N] = size(A);
delta = M/N;
ampsim_tau2 = zeros(inter_max, 1);
ampsim_MSE = zeros(inter_max, 1);

z = y;
x = zeros(N,1);
tau2 = 1/M*(norm(z,2)^2);

% figure
ampsim_tau2(1) = tau2;
ampsim_MSE(1) = norm(x_true - x,2).^2/N;
for i=1:inter_max
    x = wthresh(A'*z + x,'s',sqrt(tau2));
    
    eta_prime = 1/N * nnz(x);
    z = y - A*x + 1/delta * z * eta_prime;
%     eta_prime = 1/N * nnz(wthresh(A'*z + x,'s',sqrt(tau2)));
%     tau2 = tau2/delta * eta_prime;
    tau2 = 1/M*(norm(z,2)^2);
    
    ampsim_tau2(i+1) = tau2;
    ampsim_MSE(i+1) = norm(x_true - x,2).^2/N;
    
    if(i == 10)
        unthresholded_estimate = A'*z + x;
    end
        
%     stem(x);
%     hold on
%     stem(x_true,'r')
%     hold off
%     drawnow;
%     pause;
    fprintf('AMP MSE = %f \n', 10*log10(ampsim_MSE(i+1)));
    
end

MSE = 10*log10(norm(x_true - x,2).^2/N);
fprintf('AMP final MSE = %f ', MSE);
