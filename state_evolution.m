function [ se_tau2, se_mse ] = state_evolution(tau2, delta, rho, sigmaw2, lambda)

% state_evolution - this function outputs state evolution predictions
%
% Inputs:
% tau2      - the effective noise variance from the previous interation
% delta     - the measurement ratio (M/N)
% SigmaW    - the variance of the noise.
% rho       - the probability of a nonzero component, i.e., sparsity ratio (K/N).
%
% Outputs:
% tau2      - the effective noise variance
% mse       - the mean squared error
%
% Author:   Osman Musa
% email:    osman.musa@nt.tuwien.ac.at
% Website:  https://www.nt.tuwien.ac.at/about-us/staff/osman-musa/
% Last revision: 01-Aug-2017

if(sigmaw2 == 0)
    deniser_parameter = sqrt(tau2);
else
    deniser_parameter = sqrt(lambda+lambda*tau2);
end


% Gauss part
f1 = @(x,z) (rho.*reshape(normpdf(x(:),0,1),size(x)) ) ...
            .* reshape(normpdf(z(:),0,sqrt(tau2)),size(z)) ...
            .* (wthresh(x + z,'s',deniser_parameter) - x).^2;

% Dirac part        
f2 = @(z) (1-rho) ...
            .* reshape(normpdf(z(:),0,sqrt(tau2)),size(z)) ...
            .* (wthresh( z,'s',deniser_parameter)).^2;        
         

lim = inf;

se_mse = (integral2(f1,-lim,lim,-lim,lim) + integral(f2,-lim,lim));
se_tau2 = sigmaw2 + 1/delta * se_mse;   

%% Attempt to estimate the original denoiser parameter obtained while deriving amp
% MC      = 1e3;                    % No. of MC experiments.
% x = randn(MC,1);
% x = repmat(x, 1, MC);
% z = randn(MC,MC);
% 
% gaus_part = rho * 1/MC/MC*nnz(wthresh(x + sqrt(tau2)*z,'s',sqrt(tau2)));
% dirac_part = (1-rho)* 1/MC/MC*nnz(wthresh(sqrt(tau2)*z,'s',sqrt(tau2)));
% 
% eta_prime = gaus_part + dirac_part;
% se_tau2 = tau2/delta * eta_prime;

end
