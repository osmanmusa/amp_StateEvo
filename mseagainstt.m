function [ ampsim_tau2, ampsim_mse ] = mseagainstt(N, K, MC, M, sigmaw2, inter_max)
%% Input
% N                     Size of vector x
% K                     Sparsity level
% MC                    Number of Monte-Carlo realizations                       
% M                     Number of measurements
% sigmaw2               power of the Gaussian noise
%
% 2. Be careful with normalization of Phi when compareing different 
% algorithms. Some have different requirements for the normalization
%% Check number of inputs.


%% Find MSE for each measurement rate
tic
    
MSE_AMP = zeros(MC,1);
ampsim_tau2_single = zeros(inter_max+1, MC);
ampsim_mse_single = zeros(inter_max+1, MC);

%% Start Monte-Carlo loop % parfor
for mc=1:MC
%         disp(sprintf('-----------The %d. iteration of mc loop-----------',mc));
    %% Sensing matrix 
    A_unnormalized = sqrt(1/M)*randn(M,N);
    s = sqrt(sum(A_unnormalized.^2));
    S = diag(1./s);
    A = A_unnormalized*S; % normalized columns of A
    
%     A = sqrt(1/M)*sign(randn(M,N));

    %% Unknown data 
    ind = randperm(N);
    x = zeros(N,1);
    
    x(ind(1:K)) = randn(K,1); % create a Gaussian K-sparse source
%         x(ind(1:K)) = rand(K,1); % create a uniform 0-1 K-sparse source
%         x(ind(1:K)) = ones(K,1);
%     x = abs(x);

%     x(ind(1:K)) = sign(randn(K,1)); % create a +-1 K-sparse source
    
%     fprintf('Transmit vector power 1/M norm(Ax)^2 = %f \n', 1/M*norm(A*x)^2)
    
    noise = sqrt(sigmaw2)*randn(M,1);
    signal = A*x;
    y = signal + noise;
    
    fprintf('snr = %.3f \n',snr(signal,noise));
    
%     figure
%     plot(signal)
%     hold on
%     plot(noise,'r')
    
    %% AMP Reconstruction analysis
    if(sigmaw2 == 0)
        [x_AMP, se_tau2, se_mse] = amp0(A, y, x, inter_max);
    else
        [x_AMP, se_tau2, se_mse] = ampa(A, y, x, inter_max);
    end
    
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

