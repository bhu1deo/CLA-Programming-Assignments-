% Q2.4 a part
% We will use the tolerance as norm(A*x-b); 
% Earlier we had used relative error 

 % Increasing the System Size and Reducing the tolerance both increase the
 % number of iterations required to reach convergence

% Tolerance here is error in Axk-b; xk is the kth Jacobi Iterate 

% Error trajectory plots are shown for each case 
%%

% n =10 
disp('n is 10 here');
n = 10;
e = ones(n,1);
A = spdiags([-e 2*e -e],-1:1,n,n);
A = full(A);
b = rand(n,1);

% tol = 0.1
x0 = 0.0001*rand(n,1);      % Initial Guess 
[xjac,count_jac,ejac] = jacobi_method(A,b,x0,0.1,1000);
disp('Number of Iterations to reach tolerance n = 10,tol=0.1')
count_jac


[row,col] = size(ejac);
Xjac = 0:1:col-1;
figure()
plot(Xjac,ejac)
title('n=10, tol = 0.1')
legend('Jacobi')


% tol = 0.01
x0 = 0.0001*rand(n,1);      % Initial Guess 
[xjac,count_jac,ejac] = jacobi_method(A,b,x0,0.01,1000);
disp('Number of Iterations to reach tolerance n = 10,tol=0.01')
count_jac

[row,col] = size(ejac);
Xjac = 0:1:col-1;
figure()
plot(Xjac,ejac)
title('n=10, tol = 0.01')
legend('Jacobi')


% tol = 0.001
x0 = 0.0001*rand(n,1);      % Initial Guess 
[xjac,count_jac,ejac] = jacobi_method(A,b,x0,0.001,1000);
disp('Number of Iterations to reach tolerance n = 10,tol=0.001')
count_jac

[row,col] = size(ejac);
Xjac = 0:1:col-1;
figure()
plot(Xjac,ejac)
title('n=10, tol = 0.001')
legend('Jacobi')


% tol = 0.0001
x0 = 0.0001*rand(n,1);      % Initial Guess 
[xjac,count_jac,ejac] = jacobi_method(A,b,x0,0.0001,1000);
disp('Number of Iterations to reach tolerance n = 10,tol=0.0001')
count_jac

[row,col] = size(ejac);
Xjac = 0:1:col-1;
figure()
plot(Xjac,ejac)
title('n=10, tol = 0.0001')
legend('Jacobi')


% tol = 0.00001
x0 = 0.0001*rand(n,1);      % Initial Guess 
[xjac,count_jac,ejac] = jacobi_method(A,b,x0,0.00001,1000);
disp('Number of Iterations to reach tolerance n = 10,tol=0.00001')
count_jac

[row,col] = size(ejac);
Xjac = 0:1:col-1;
figure()
plot(Xjac,ejac)
title('n=10, tol = 0.00001')
legend('Jacobi')

%%
disp('n is 50 here'); 
n = 50;
e = ones(n,1);
A = spdiags([-e 2*e -e],-1:1,n,n);
A = full(A);
b = rand(n,1);
% tol = 0.1
x0 = 0.0001*rand(n,1);      % Initial Guess 
[xjac,count_jac,ejac] = jacobi_method(A,b,x0,0.1,10000);
disp('Number of Iterations to reach tolerance n = 50,tol=0.1')
count_jac


[row,col] = size(ejac);
Xjac = 0:1:col-1;
figure()
plot(Xjac,ejac)
title('n=50, tol = 0.1')
legend('Jacobi')


% tol = 0.01
x0 = 0.0001*rand(n,1);      % Initial Guess 
[xjac,count_jac,ejac] = jacobi_method(A,b,x0,0.01,10000);
disp('Number of Iterations to reach tolerance n = 50,tol=0.01')
count_jac

[row,col] = size(ejac);
Xjac = 0:1:col-1;
figure()
plot(Xjac,ejac)
title('n=50, tol = 0.01')
legend('Jacobi')


% tol = 0.001
x0 = 0.0001*rand(n,1);      % Initial Guess 
[xjac,count_jac,ejac] = jacobi_method(A,b,x0,0.001,10000);
disp('Number of Iterations to reach tolerance n = 50,tol=0.001')
count_jac

[row,col] = size(ejac);
Xjac = 0:1:col-1;
figure()
plot(Xjac,ejac)
title('n=50, tol = 0.001')
legend('Jacobi')


% tol = 0.0001
x0 = 0.0001*rand(n,1);      % Initial Guess 
[xjac,count_jac,ejac] = jacobi_method(A,b,x0,0.0001,10000);
disp('Number of Iterations to reach tolerance n = 50,tol=0.0001')
count_jac

[row,col] = size(ejac);
Xjac = 0:1:col-1;
figure()
plot(Xjac,ejac)
title('n=50, tol = 0.0001')
legend('Jacobi')


% tol = 0.00001
x0 = 0.0001*rand(n,1);      % Initial Guess 
[xjac,count_jac,ejac] = jacobi_method(A,b,x0,0.00001,10000);
disp('Number of Iterations to reach tolerance n = 50,tol=0.00001')
count_jac

[row,col] = size(ejac);
Xjac = 0:1:col-1;
figure()
plot(Xjac,ejac)
title('n=50, tol = 0.00001')
legend('Jacobi')


%%
disp('n is 100 here'); 
n = 100;
e = ones(n,1);
A = spdiags([-e 2*e -e],-1:1,n,n);
A = full(A);
b = rand(n,1);
% tol = 0.1
x0 = 0.0001*rand(n,1);      % Initial Guess 
[xjac,count_jac,ejac] = jacobi_method(A,b,x0,0.1,100000);
disp('Number of Iterations to reach tolerance n = 100,tol=0.1')
count_jac


[row,col] = size(ejac);
Xjac = 0:1:col-1;
figure()
plot(Xjac,ejac)
title('n=100, tol = 0.1')
legend('Jacobi')


% tol = 0.01
x0 = 0.0001*rand(n,1);      % Initial Guess 
[xjac,count_jac,ejac] = jacobi_method(A,b,x0,0.01,100000);
disp('Number of Iterations to reach tolerance n = 100,tol=0.01')
count_jac

[row,col] = size(ejac);
Xjac = 0:1:col-1;
figure()
plot(Xjac,ejac)
title('n=100, tol = 0.01')
legend('Jacobi')


% tol = 0.001
x0 = 0.0001*rand(n,1);      % Initial Guess 
[xjac,count_jac,ejac] = jacobi_method(A,b,x0,0.001,100000);
disp('Number of Iterations to reach tolerance n = 100,tol=0.001')
count_jac

[row,col] = size(ejac);
Xjac = 0:1:col-1;
figure()
plot(Xjac,ejac)
title('n=100, tol = 0.001')
legend('Jacobi')


% tol = 0.0001
x0 = 0.0001*rand(n,1);      % Initial Guess 
[xjac,count_jac,ejac] = jacobi_method(A,b,x0,0.0001,100000);
disp('Number of Iterations to reach tolerance n = 100,tol=0.0001')
count_jac

[row,col] = size(ejac);
Xjac = 0:1:col-1;
figure()
plot(Xjac,ejac)
title('n=100, tol = 0.0001')
legend('Jacobi')


% tol = 0.00001
x0 = 0.0001*rand(n,1);      % Initial Guess 
[xjac,count_jac,ejac] = jacobi_method(A,b,x0,0.00001,100000);
disp('Number of Iterations to reach tolerance n = 100,tol=0.00001')
count_jac

[row,col] = size(ejac);
Xjac = 0:1:col-1;
figure()
plot(Xjac,ejac)
title('n=100, tol = 0.00001')
legend('Jacobi')
% 




%%





% Jacobi Method with tolerance as stated above 
% If method doesn't converge stop at some MAX iterations 
function [xjac,count_jac,ejac] = jacobi_method(A,b,x0,precision,maxit)
format
% first check whether convergence possible from
% arbitrary initialization
 

L = zeros(3);
U = zeros(3);
D = zeros(3);

L(2:3,1) = A(2:3,1);
L(3,2) = A(3,2);
U(1:2,3) = A(1:2,3);
U(1,2) = A(1,2);
D(1,1) = A(1,1);
D(2,2) = A(2,2);
D(3,3) = A(3,3);

PJ = -inv(D)*(L + U);       % P JAcobi



disp('Spectral radius of P Jacobi')
lambdamaxj = max(abs(eig(PJ)))


if(lambdamaxj<1)
    disp('Convergence for Jacobi possible from arbitrary initialization')
else
    disp('Convergence for Jacobi not possible from arbitrary initialization')
end

% Run both algorithms either till precision is met or max no. of iterations 
% Assumed that an nxn invertible matrix is passed here and diags are
% nonzero
C = A;
n1 = size(A);
n = n1(1);                               % Rows = Columns here

 
% First the Jacobi Method, here relative error is used as precision
% norm(A*x-b)>epsilon;  % This also we can use alternatively

x = x0;
epsilon = precision;
count_jac = 0;
ejac = 0;
while(norm(A*x-b)>epsilon & count_jac<maxit)
    X = repmat(x,1,n);           % replicate the current iter
    X = X - diag(diag(X));       % i,i to be made 0
    xtemp = diag(A*X);           % pick only diagonal values
    xtemp = -xtemp./diag(A);     % divide by A(i,i)
    ejac(end+1) = norm(xtemp + b./diag(A)-x)/(norm(xtemp + b./diag(A))+0.0000001);
    x = xtemp + b./diag(A);      % add the b term divide by A(i,i) 
    count_jac=count_jac+1;
end
xjac = x;

end





