% Here we do Householder reflection to compute the QR decomposition 
% Index in MATLAB starts from 1
  

function [Q,R] = hr(A)
format 
n1 = size(A);
n = n1(1);                               % Rows 
m = n1(2);                  % The number of columns in A


B = A;                        % The matrix on which Householder process is applied 
R = zeros(n,m);              % For storing the R matrix 

for k = 1:m
   x = B(:,1);                  % The vector under consideration 
   x1 = size(x);
   xr = x1(1);                 % number of rows 
   y = zeros(xr,1);
   y(1) = norm(x);
   u = x-y;                     % The per vector 
   if(norm(u)~=0)                
       u = u/norm(u);
   end
   
                 
   
   Q = eye(xr)-2*u*transpose(u);
   B = Q*B;
  
   
   if(k==1)
       Qf=Q;
   else
       
       Qt = eye(n);                   % That Submatrix thing 
       Qt(n-xr+1:end,n-xr+1:end) = Q;
       Qf = Qt*Qf;
   end
       
   R(k:n,k) = B(:,1);               % Add the column 
   if(k<=n)
       R(k,k:m) = B(1,:);                % Add the row 
   end
   if(k==m)
       break
   end
   B = B(2:end,2:end);             % Now operate on this 
end

Q = transpose(Qf);
disp("Q")
disp(Q)
disp("R")
disp(R)
disp("Q*R")
disp(Q*R)
disp("A")
disp(A)
disp("QTQ")
disp(transpose(Q)*Q)
disp('Frobenius NORM of the error matrix')
disp(norm(A-Q*R))
end
