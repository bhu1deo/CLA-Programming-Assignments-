% Gram Schmidt process to get QR decomposition 
% This is the classical Gram Schmidt process
% Index in MATLAB starts from 1
% Here we will do QR using Gram Schmidt for the LI case only 
% This code works for tall as well as fat matrices 

function [Q,R] = gs(A)          % Return Q and R for a given A
format 
n = size(A);
m = n(2);                  % The number of columns in A
Q = zeros(n(1),n(2));         % Q init
R = zeros(n(2),n(2));                   % R init
 
for k=1:m
    if(k==1)
        Q(:,1) = A(:,1)/norm(A(:,1));          % 1st vector 
        R(1,1) = norm(A(:,1));
    else
        v = A(:,k);                  % The current vector
        r = transpose(v)*Q;           % The projections 
        
        q = v - Q*transpose(r);       % Component ortho to all prev proj.
        r(k) = norm(q);               % R[k,k]
        Q(:,k) = q/norm(q);           % Orthonormality
        R(:,k) = transpose(r);         % The kth col. of R
    end

end

disp("The Original Matrix")
disp(A)
disp("Q")
disp(Q)
disp("R")
disp(R)
disp("Q*R")
disp(Q*R)
disp("QTQ")
disp(transpose(Q)*Q)
disp('Frobenius NORM of the error matrix')
disp(norm(A-Q*R))
end



