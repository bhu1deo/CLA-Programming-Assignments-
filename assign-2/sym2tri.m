%Write a program (sym2tri.m) to convert a symmetric matrix to tridiagonal
%matrix using a
%similarity transformation. % Input is a symmetric matrix. Output is a
%tridiagonal matrix. (Special case of Hessenberg Matrix)
% Input is assumed to be a square symmetric matrix 

function [tridiag] = sym2tri(A)
if(A~=transpose(A))
    disp('Matrix is not symmetric, tridiag form not possible')
    return 
end

% Matrix is Symmetric case below:
n = size(A);
n = n(1);
% We don't need Q so we don't return it 
% We will not do a signed implementation rather check if the column is as
% desired then do nothing NORM is preserved under reflection
B = A;            % Keep operating on this instead 
for k=1:n-2 
    x = B(k+1:end,k);        % The current column
    if(all(x(2:end,1))==0)     % Nothing to do 
        continue
    end
    temp = zeros(n-k,1);
    
    temp(1) = norm(x);
    u = x-temp;
    u = u/norm(u);     % This is non zero
    % First calculate the ortho matrix and then embed it into identity
    % matrix 
    Qstar = eye(n-k)-2*u*transpose(u);
    % Now embed it 
    Q = eye(n);
    Q(k+1:end,k+1:end) = Qstar;
    B = Q*B*transpose(Q);% We could explicitly avoid this but is done for completeness
end
tridiag = B;             % Final tridiagonal form 
end