% Q 2.4 b part CLA 

% Our best performing direct method is just the QR with rotation and then
% solving Ry = QTb; when initially it was QRx = b; Then we have to compare
% with Jacobi methods with given tolerance and different values of n =
% 10;50;100. Just time required should be compared 

% In general best performing direct method gave better precision and lesser
% computational time was required 

% Better to run the codes sectionwise 

%%
disp('n is 10');
n=10;
e = ones(n,1);
A = spdiags([-e 2*e -e],-1:1,n,n);
A = full(A);
b = rand(n,1);

precision = 0.001;       % for jacobi method 
x0 = 0.0001*rand(n,1);    % for jacobi method 

[tjac,xjac] = solve_jacobi(A,b,x0,precision); 
disp('Jacobi method precision') 
disp(norm(A*xjac-b));
disp('Jacobi Method Time Required')
disp(tjac)
[tdirect,xdirect] = solve_QR(A,b);              % Direct Method 
disp('QR PROP method precision') 
disp(norm(A*xdirect-b));
disp('QR PROP Method Time Required')
disp(tdirect);
%%
disp('n is 50');
n=50;
e = ones(n,1);
A = spdiags([-e 2*e -e],-1:1,n,n);
A = full(A);
b = rand(n,1);

precision = 0.001;       % for jacobi method 
x0 = 0.0001*rand(n,1);    % for jacobi method 

[tjac,xjac] = solve_jacobi(A,b,x0,precision); 
disp('Jacobi method precision') 
disp(norm(A*xjac-b));
disp('Jacobi Method Time Required')
disp(tjac)
[tdirect,xdirect] = solve_QR(A,b);              % Direct Method 
disp('QR PROP method precision') 
disp(norm(A*xdirect-b));
disp('QR PROP Method Time Required')
disp(tdirect);
%%
disp('n is 100');
n=100;
e = ones(n,1);
A = spdiags([-e 2*e -e],-1:1,n,n);
A = full(A);
b = rand(n,1);

precision = 0.001;       % for jacobi method 
x0 = 0.0001*rand(n,1);    % for jacobi method 

[tjac,xjac] = solve_jacobi(A,b,x0,precision); 
disp('Jacobi method precision') 
disp(norm(A*xjac-b));
disp('Jacobi Method Time Required')
disp(tjac)
[tdirect,xdirect] = solve_QR(A,b);              % Direct Method 
disp('QR PROP method precision') 
disp(norm(A*xdirect-b));
disp('QR PROP Method Time Required')
disp(tdirect);
%%
function [time_req,x_jac] = solve_jacobi(A,b,x0,precision)
% return the time of computation to reach the required precision 
tic;                % Start the time computation 
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

PJ = -inv(D)*(L + U);       % P Jacobi

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
ejac = 0;
while(norm(A*x-b)>epsilon)
    X = repmat(x,1,n);           % replicate the current iter
    X = X - diag(diag(X));       % i,i to be made 0
    xtemp = diag(A*X);           % pick only diagonal values
    xtemp = -xtemp./diag(A);     % divide by A(i,i)
    ejac(end+1) = norm(xtemp + b./diag(A)-x)/(norm(xtemp + b./diag(A))+0.0000001);
    x = xtemp + b./diag(A);      % add the b term divide by A(i,i) 
end

x_jac = x;

time_req = toc;         % Time required to finish the computation 
end

function [time_req,x_direct] = solve_QR(A,b)

% return the time required to solve using direct QR method 
% This solution would be exact limited to computation errors 
tic;                   % Start time computation 
% First do the QR decomposition; QRx = b; Rx = QTb 
n = size(A);
n = n(1);                       % Square matrix that's why
[Q,R] = prop(A);

% Then solve for Rx = b' = QTb
b_prime = transpose(Q)*b;          % changed b 
x_vec = zeros(n,1);
for k=1:n
   index = n-k+1;
   if(index==n)
       x_vec(index) = b_prime(index)/R(index,index);
   else
       x_vec(index) = (b_prime(index)-R(index,index+1:end)*x_vec(index+1:end,1))/R(index,index);
   end
end
x_direct = x_vec; 
time_req = toc; 

end

