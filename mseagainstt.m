function [ ampsim_tau2, ampsim_mse ] = mseagainstt(N, K, MC, M, sigmaw2, lambda, inter_max)

% mseagainstt - performs monte-carlo simulation to get averaged effective
% variance and mse accross iterations
%
% Inputs:
% N         - size of vector x
% K         - sparsity level
% MC        - number of Monte-Carlo realizations                       
% M         - number of measurements
% sigmaw2   - variance of the Gaussian noise
% lambda    - denoiser parameter in the noisy case
% inter_max - maximum number of amp iterations
%
% Outputs:
% tau2      - vector of averaged effective variances accross interations
% mse       - vector of averaged mean squared errors accross interations
%
% Author:   Osman Musa
% email:    osman.musa@nt.tuwien.ac.at
% Website:  https://www.nt.tuwien.ac.at/about-us/staff/osman-musa/
% Last revision: 01-Aug-2017

%% Find MSE for each measurement rate
tic
    
ampsim_tau2_single = zeros(inter_max, MC);
ampsim_mse_single = zeros(inter_max, MC);

%% Start Monte-Carlo loop 
for mc=1:MC % parfor
    %% Sensing matrix 
    A_unnormalized = sqrt(1/M)*randn(M,N);
    s = sqrt(sum(A_unnormalized.^2));
    S = diag(1./s);
    A = A_unnormalized*S; % normalized columns of A

    %% Unknown sparse vector 
    x = zeros(N,1);
    ind = randperm(N);
    
    x(ind(1:K)) = randn(K,1);           % create a Gaussian K-sparse source
%     x(ind(1:K)) = rand(K,1);          % create a uniform 0-1 K-sparse source
%     x(ind(1:K)) = ones(K,1);          % create a +1 K-sparse source
%     x(ind(1:K)) = sign(randn(K,1));   % create a +-1 K-sparse source
    %% Measurement process
    
    noise = sqrt(sigmaw2)*randn(M,1);
    signal = A*x;
    y = signal + noise;

    fprintf('%d. iteration: snr = %.3f \n', mc, snr(signal,noise)); % optional printing

    %% AMP Reconstruction
    [~, se_tau2, se_mse] = amp(A, y, lambda, x, (sigmaw2 == 0), inter_max);
    
    ampsim_tau2_single(:,mc) = se_tau2;
    ampsim_mse_single(:,mc) = se_mse;
end % end-for Monte-Carlo loop


%% Average over tau2 and mse and save results
ampsim_tau2 = mean(ampsim_tau2_single,2);
ampsim_mse = mean(ampsim_mse_single,2);

save(sprintf('SE K=%d',K),'ampsim_tau2','ampsim_mse');
toc

end

