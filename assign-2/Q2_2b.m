% Q 2.2b test the CorrNRV function 

mu1 = [0;0];
sigma1 = [0.025,0.0075;0.0075,0.007];

mu2 = [0;0;0];
sigma2 = [0.025,0.0075,0.00175;0.0075,0.007,0.00135;0.00175,0.00135,0.00043];

[Y1,x1] = corrNRV(1000,2,mu1,sigma1);

[Y2,x2] = corrNRV(1000,3,mu2,sigma2);

% x1 contains 1000 2D uncorrelated vectors Y1 contains 1000 correlated
% vectors with given sigma1

% x2 contains 1000 3D uncorrelated vectors Y2 contains 1000 correlated
% vectors with given sigma2
