
load('variables')
histogram(reshape(AMP_unthresholded_residuals,MC*K,1), 100,'Normalization','pdf')
y = fitdist(reshape(AMP_unthresholded_residuals,MC*K,1), 'Normal');
hold on
x = -1:0.01:1;
f = normpdf(x, y.mu, y.sigma);
plot(x,f,'r')