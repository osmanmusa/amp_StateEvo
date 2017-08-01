%% Simulation parameters
vector_of_K = [5 50 500];
MC_simulation = 2;
inter_max = 20;
snr_proxy = 0.0;

% allocating zero vectors
se_tau2 = zeros(inter_max,1);
se_mse = zeros(inter_max,1);

figure

%% Plotting graphs for each K
for index = 1:length(vector_of_K)
    
    %% Problem parameters
    K = vector_of_K(index);
    N = 4000;
    M = 800;
    delta = M/N;
    rho = K/N;
    sigmaw2 = snr_proxy*K/N; %K/N is the power of x
    lambda = 0.1;
    
    %% Numerical simulation
    % starts a new set of simulatios
    [ampsim_tau2, ampsim_mse] = mseagainstt(N, K, MC_simulation, M, sigmaw2, lambda, inter_max); 
    % loads the simulation results from a file, e.g. "SE K=50.mat"
%     load(sprintf('SE K=%d',K));
    
    %% SE prediction
    se_mse(1) = K/N;
    se_tau2(1) = sigmaw2 + 1/delta*se_mse(1);
    for i=2:inter_max
        [ se_tau2(i), se_mse(i) ] = state_evolution(se_tau2(i-1), delta, rho, sigmaw2, lambda);
    end
    
    %% plots in dB
    subplot(length(vector_of_K),2,2*index-1)
    plot([0:inter_max-1], 10*log10(ampsim_tau2),'m-');
    hold on
    plot([0:inter_max-1], 10*log10(ampsim_mse),'c-');
    plot([0:inter_max-1], 10*log10(se_tau2),'r*');
    plot([0:inter_max-1], 10*log10(se_mse),'bx');
    
    ylabel('MSE [dB]')
    xlabel('iteration #')
    title(sprintf('K=%d, rho=%.3f',K,K/N))
    %% plots in magnitude    
    subplot(length(vector_of_K),2,2*index)  
    plot([0:inter_max-1], ampsim_tau2,'m-');
    hold on
    plot([0:inter_max-1], ampsim_mse,'c-');
    plot([0:inter_max-1], se_tau2,'r*');
    plot([0:inter_max-1], se_mse,'bx');

    legend('Simulation tau2', 'Simulation MSE', 'SE tau2', 'SE MSE')
    ylabel('MSE')
    xlabel('iteration #')
    title(sprintf('K=%d, rho=%.3f',K,K/N))
end