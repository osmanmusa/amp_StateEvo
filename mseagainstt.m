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
    
MSE_AMP = zeros(MC,1);
ampsim_tau2_single = zeros(inter_max, MC);
ampsim_mse_single = zeros(inter_max, MC);

%% Start Monte-Carlo loop % parfor
for mc=1:MC
%         disp(sprintf('-----------The %d. iteration of mc loop-----------',mc));

    %% Sensing matrix 
    A_unnormalized = sqrt(1/M)*randn(M,N);
    s = sqrt(sum(A_unnormalized.^2));
    S = diag(1./s);
    A = A_unnormalized*S; % normalized columns of A

    %% Unknown data 
    ind = randperm(N);
    x = zeros(N,1);
    
    x(ind(1:K)) = randn(K,1); % create a Gaussian K-sparse source
%         x(ind(1:K)) = rand(K,1); % create a uniform 0-1 K-sparse source
%         x(ind(1:K)) = ones(K,1);
%     x = abs(x);

%     x(ind(1:K)) = sign(randn(K,1)); % create a +-1 K-sparse source
    %% Measurement process
    noise = sqrt(sigmaw2)*randn(M,1);
    signal = A*x;
    y = signal + noise;
    
    %     fprintf('Signal power 1/M norm(Ax)^2 = %f \n', 1/M*norm(signal)^2) 
    fprintf('%d. iteration: snr = %.3f \n', mc, snr(signal,noise));
    
%     figure
%     plot(signal)
%     hold on
%     plot(noise,'r')
    
    %% AMP Reconstruction
    amp0 = (sigmaw2 == 0);
    [x_AMP, se_tau2, se_mse] = amp(A, y, lambda, x, amp0, inter_max);
    
    MSE_AMP(mc) = 1/N*norm(x - x_AMP,2).^2;
    ampsim_tau2_single(:,mc) = se_tau2;
    ampsim_mse_single(:,mc) = se_mse;
end % end-for Monte-Carlo loop


%% Normalized MSE (in dB)
Results.MSE_AMP = 10*log10(mean(MSE_AMP));
fprintf('MSE = %f ', Results.MSE_AMP);

ampsim_tau2 = mean(ampsim_tau2_single,2);
ampsim_mse = mean(ampsim_mse_single,2);


save(sprintf('SE K=%d',K),'ampsim_tau2','ampsim_mse');

sum_time = toc
end

