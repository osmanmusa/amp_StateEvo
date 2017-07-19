function [ x, nit, unthresholded_estimate ] = ist(Phi, y, lambda, x_true)
% function [x , nit ] = IST1(Phi,y,lambda)
%     
% This iterative soft thresholding with lambda-regularization.
%
%  INPUTS :
%  Phi    :  measurement matrix m x n
%  y      :  measurement vector m x 1
%  lambda : regularization parameter, 1 is a often good choice
%
%  
%  OUTPUTS:
%  x      : recovered signal vector
%  nit    : number of iterations required 
%
%
%  
% NG, 15 May 2014 
%
% -----------------------------------------------------------------------

szP = size(Phi);

% data samples
n = szP(2);

% measurements
m = szP(1);
sqm = sqrt(m);

% Set max number of iteations to max  factor times signal dimension n
% MaxIT = 40 * n;
MaxIT = 40;
% ---------------------
% Initialisation:
z = y;
x = zeros(n,1);
xold = zeros(n,1);

PhiT = Phi';

% --------------------

% Scaling of the matrix: this is still needed, even though the Phi-matrix
% has indeed normalised columns
mu = 2/((norm(Phi))^2);

figure
j= 0;
while ( (j<MaxIT) )
    
    j = j + 1;
                     
    w =  x +  mu .* PhiT * z;
                      
    % Soft thresholding:
    x = sign(w)  .* max( abs(w) - lambda, 0);
                    
    z = y - Phi * x;   
    
    if(j == 10)
        unthresholded_estimate = Phi'*z + x;
    end
     
    % Stopping criterion
    normx = norm(x);
    normdx = norm(x-xold);
    
    if(  normdx < 1E-4 * normx )
        break;
    end
    
    xold = x;     
    
            
    stem(x);
    hold on
    stem(x_true,'r')
    hold off
    drawnow;
    pause;
    
end

x    = x(:);
nit  = j;

MSE = 10*log10(norm(x_true - x,2).^2);
message = sprintf('IST MSE = %f ', MSE);
disp(message);

end





