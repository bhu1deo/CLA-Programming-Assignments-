% corrNRV.m generate correlated samples from uncorrelated samples 


function [Y,x] = corrNRV(M,n,mu,sigma)      % M,n mu and sigma are taken as arguments 

x = randn(n,M);  % M random vectors each of n dimension are generated from N(0,I);
L = cholprog(sigma);    % The dimension of the sigma matrix should be nxn    
% This is because the dimension of each random vector is nx1, dimension of
% L is also nxn hence mu is vector of dimension nx1 
mu_new = repmat(mu,1,M);      % Note that mu is a nx1 vector 
Y = L*x + mu_new;      % This is the required correlated data matrix        

end