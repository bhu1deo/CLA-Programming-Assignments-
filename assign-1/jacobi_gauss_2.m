% Q2.3.2
% We will use the same function jacobi_gauss as earlier 
%%
% First case 
r = zeros(10,1);
r(1) = 2;
r(2) = -1;
A = toeplitz(r);      % A is strictly diagonally dominant
% We thus don't need to compute the spectral radius of P,Jacobi and GS will
% converge 
b = zeros(10,1);
b(1)=11;
x0 = zeros(10,1);
maxit = 1000; 
precision = 0.0001;              % precision 
[xjac,xgs,count_jac,count_gs,ejac,egs] = jacobi_gauss(A,b,x0,precision,maxit);
% A is strictly diagonally dominant
% We thus don't need to compute the spectral radius of P,Jacobi and GS will
% converge

%%

% Second Part 
r = zeros(20,1);
r(1) = 2;
r(2) = -1;
A = toeplitz(r);      % A is strictly diagonally dominant
% We thus don't need to compute the spectral radius of P,Jacobi and GS will
% converge 
b = zeros(20,1);
b(1)=21;
x0 = zeros(20,1);
maxit = 1000; 
precision = 0.0001;              % precision 
[xjac,xgs,count_jac,count_gs,ejac,egs] = jacobi_gauss(A,b,x0,precision,maxit);

% A is strictly diagonally dominant
% We thus don't need to compute the spectral radius of P,Jacobi and GS will
% converge

%%

function [xjac,xgs,count_jac,count_gs,ejac,egs] = jacobi_gauss(A,b,x0,precision,maxit)
format
% For both the methods first check whether convergence possible from
% arbitrary initialization

xstar = inv(A)*b;       % Only used to compute the relative error
if(x0==xstar)
    disp('Nothing to do, you have entered the optimal solution')
    xjac = xstar;
    xgs = xstar;
    count_jac = 0;
    count_gs = 0;
    ejac = 0;
    egs = 0;
    return;
end

init_error = norm(x0-xstar);          % Used to compute relative error 

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

PJ = -inv(D)*(L + U);       % P Jacobi
PGS = -inv(L+D)*U;                % P Gauss Seidel 

disp('Spectral radius of P Gauss Seidel')
lambdamax_gs = max(abs(eig(PGS)))
disp('Spectral radius of P Jacobi')
lambdamaxj = max(abs(eig(PJ)))

if(lambdamax_gs<1)
    disp('Convergence for Gauss Seidel possible from arbitrary initialization')
else
    disp('Convergence for Gauss Seidel not possible from arbitrary initialization')
end

if(lambdamaxj<1)
    disp('Convergence for Jacobi possible from arbitrary initialization')
else
    disp('Convergence for Jacobi not possible from arbitrary initialization')
end

% Run both algorithms either till precsion is met or max no. of iterations 
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
while(norm(x-xstar)/init_error>epsilon & count_jac<maxit)
    X = repmat(x,1,n);           % replicate the current iter
    X = X - diag(diag(X));       % i,i to be made 0
    xtemp = diag(A*X);           % pick only diagonal values
    xtemp = -xtemp./diag(A);     % divide by A(i,i)
    ejac(end+1) = norm(xtemp + b./diag(A)-x)/(norm(xtemp + b./diag(A))+0.0000001);
    x = xtemp + b./diag(A);      % add the b term divide by A(i,i) 
    count_jac=count_jac+1;
end

xjac = x;

% Now the Gauss Seidel Method, this is a bit difficult to vectorize 
% We will do normally double for loops 

x = x0;            % The starting point 
count_gs = 0;
egs = 0;
while(norm(x-xstar)/init_error>epsilon & count_gs<maxit )
    count_gs=count_gs+1;
    xprev = x;
    for j=1:n
       if(j==1)          % only the last iterate
           x(1,:) = -A(1,2:n)*x(2:n,1)+b(1,1);
           x(1,:) = x(1,:)/A(1,1);
           
       elseif(j==n)         % completely this iterate
           x(n,:) = -A(n,1:n-1)*x(1:n-1,1)+b(n,1);
           x(n,:) = x(n,:)/A(n,n);
        
       else                % intermediate
           x(j,:) = -A(j,1:j-1)*x(1:j-1,1)-A(j,j+1:n)*x(j+1:n,1)+b(j,1);
           x(j,:) = x(j,:)/A(j,j);
           
       end
       
    end
    
    egs(end+1) = norm(x-xprev)/(norm(x)+0.000000001);
    
end
ejac = ejac(2:end);
egs = egs(2:end);
% Now Plot it 
% First the Jacobi Method error trajectory 
[row,col] = size(ejac);
Xjac = 0:1:col-1;
[row,col] = size(egs);
Xgs = 0:1:col-1;
plot(Xjac,ejac,Xgs,egs)
title('Rate of Convergence')
legend('Jacobi','Gauess Seidel')
% Now the Gauss Seidel method error trajectory
xgs = x;                    % The gauss seidel iterate 


end

