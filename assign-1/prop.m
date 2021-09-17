% Doing QR decomposition using rotation matrices 
% A rotator matrix is the perturbed identity matrix having the ith and the 
% jth rows changed, it affects only the ith and the jth entries of the
% vector, rest entries remain unchanged 
%  here I do QR based on rot. mat.   

function [Q,R] = prop(A)
format 
n1 = size(A);
n = n1(1);                               % Rows 
m = n1(2);                  % The number of columns in A
B = A;                       % Keep the original matrix intact 
Q = eye(n);                    % Initialize to identity matrix 
for k = 1:m
    % process the current column
  
    for j = (k+1):n
        x = B(:,k);               % The kth column 
        Qt = eye(n); 
        if(x(k)==0 && x(j)==0)
            c = 1;                            % cos(theta)
            s = 0;                            % sin(theta)
        else
            beta = max(abs(x(k)),abs(x(j)));        % clearly beta won't be 0
            xcap_k = x(k)/beta;
            xcap_j = x(j)/beta;
            new = (xcap_k^2 + xcap_j^2)^0.5;
            c = xcap_k/new;
            s = xcap_j/new;                        % cos and sin
        Qt(k,k) = c;
        Qt(j,j) = c;
        Qt(j,k) = -1*s;
        Qt(k,j) = s; 
        Q = Qt*Q;           % To Compute the final Q matrix 
        B = Qt*B;              % One index pair rotated
        end
        
    end
   
end
R = B;
disp('A')
A
disp('Q')
Q = transpose(Q)
disp('R')
R
disp('Error')
norm(A-Q*R)
disp('QTQ')
transpose(Q)*Q
end











