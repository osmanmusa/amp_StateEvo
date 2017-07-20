function [ se_tau2, se_mse ] = func_SE_MC(tau2, delta, rho, sigmaw2)
% SigmaX is the variance of the nonzero components in the signal.
% SigmaW is the variance of the effective noise.
% rho is the probability of a nonzero component, i.e., sparsity ratio (K/N).
  
if(sigmaw2 == 0)
    deniser_parameter = sqrt(tau2);

else
    lambda = 0.5;
    deniser_parameter = sqrt(lambda+lambda*tau2);
end


% Gauss part
f1 = @(x,z) (rho.*reshape(normpdf(x(:),0,1),size(x)) ) ...
            .* reshape(normpdf(z(:),0,deniser_parameter),size(z)) ...
            .* (wthresh(x + z,'s',deniser_parameter) - x).^2;

% Dirac part        
f2 = @(z) (1-rho) ...
            .* reshape(normpdf(z(:),0,deniser_parameter),size(z)) ...
            .* (wthresh( z,'s',deniser_parameter)).^2;        


lim = inf;

se_mse = (integral2(f1,-lim,lim,-lim,lim) + integral(f2,-lim,lim));
se_tau2 = sigmaw2 + 1/delta * se_mse;   


end
