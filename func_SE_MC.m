function [ se_tau2, se_mse ] = func_SE_MC(tau2, delta, rho, sigma2)
% SigmaX is the variance of the nonzero components in the signal.
% SigmaW is the variance of the effective noise.
% rho is the probability of a nonzero component, i.e., sparsity ratio (K/N).
  
MC      = 1e6;                    % No. of MC experiments.

ind = randperm(MC);
x = zeros(MC,1);
x(ind(1:round(MC*rho,0))) = randn(round(MC*rho,0),1);
z = randn(MC,1);

if(sigma2 == 0)
    se_mse = 1/MC*sum((wthresh(x + sqrt(tau2)*z,'s',sqrt(tau2)) - x).^2);
    se_tau2 = sigma2 + 1/delta * se_mse;
else
    lambda = 0.5;
    se_mse = 1/MC*sum((wthresh(x + sqrt(tau2)*z,'s',sqrt(lambda+lambda*tau2)) - x).^2);
    se_tau2 = sigma2 + 1/delta * se_mse;
end