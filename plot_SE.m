function [] = plot_SE()

N = 4000;
M = 2000;
delta = M/N;

Ks = [50 500];
MC_simulation = 10;
inter_max = 10;

se_tau2 = zeros(inter_max+1,1);
se_mse = zeros(inter_max+1,1);

figure
for index = 1:2
    
    K = Ks(index);
    rho = K/N;
    sigmaw2 = 0.0*1*K/N;
    
    %% Numerical simulation
    % starts a new set of simulatios
    [ampsim_tau2, ampsim_mse] = mseagainstt(N, K, MC_simulation, M, sigmaw2, inter_max); 
    % loads the simulation results from a file, e.g. "SE K=50.mat"
%     load(sprintf('SE K=%d',K));
    
    %% SE prediction
    se_mse(1) = K/N;
    se_tau2(1) = sigmaw2 + 1/delta*se_mse(1);
    for i=2:inter_max+1
        [ se_tau2(i), se_mse(i) ] = func_SE_MC(se_tau2(i-1), delta, rho, sigmaw2);
    end
 
    plot_helper(K, N, ampsim_tau2, ampsim_mse, se_tau2, se_mse, index, 1, inter_max); %plots in dB
    plot_helper(K, N, ampsim_tau2, ampsim_mse, se_tau2, se_mse, index, 2, inter_max); %plots in magnitude
    
end
end

function [] = plot_helper(K, N, ampsim_tau2, ampsim_mse, se_tau2, se_mse, row_index, col_index, inter_max)

    if(col_index == 1)      
        ylabel('MSE [dB]')
        subplot(2,2,2*row_index-1)
        plot([0:inter_max], 10*log10(ampsim_tau2),'m-');
        hold on
        plot([0:inter_max], 10*log10(ampsim_mse),'c-');
        plot([0:inter_max], 10*log10(se_tau2),'r*');
        plot([0:inter_max], 10*log10(se_mse),'bx');
    else
        ylabel('MSE')
        subplot(2,2,2*row_index)  
        plot([0:inter_max], ampsim_tau2,'m-');
        hold on
        plot([0:inter_max], ampsim_mse,'c-');
        plot([0:inter_max], se_tau2,'r*');
        plot([0:inter_max], se_mse,'bx');
        
        legend('Simulation tau2', 'Simulation MSE', 'SE tau2', 'SE MSE')
    end
    
    xlabel('iteration #')
    title(sprintf('K=%d, rho=%.3f',K,K/N))
    
end