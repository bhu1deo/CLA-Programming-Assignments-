% LU Decomposition without partial pivoting, computation using real numbers only 
% It's specifically given in the question that we should do for 
% square matrices, also LU doesn't exist condition needs to be 
% checked 
% All Leading principle submatrices should be nonsingular for LU to exist 

function [L,U] = cmplx(A)
format 

if(det(A)==0)           % Matrix is nonsingular 
    disp('LU not possible, matrix not invertible')
    return;
end

n1 = size(A);
n = n1(1);                               % Rows = Columns here
% Check the leading principle submatrix nonsingularity condition
for i=1:n
    if(det(A(1:i,1:i))==0)
        disp('LU not possible, do LU with pivoting')
        return;
    end
end


L = eye(n);
U = zeros(n);
U(1,:) = A(1,:);               % First row of U is same as first row of A
x = A(:,1)/U(1,1);            % The first column of L
L(:,1) = x;
% Now determine the kth row of U and the kth column of L
% interweavingly 
for k = 2:n
    
    U(k,k:n) = A(k,k:n) - L(k,1:k-1)*U(1:k-1,k:n);
    L(k+1:n,k) = (1/U(k,k))*(A(k+1:n,k)-L(k+1:n,1:k-1)*U(1:k-1,k));
end
disp('A is as follows:')
A
disp('L*U:')
L*U
disp('Frobenius Norm of the error')
norm(A-L*U)
disp('L')
L
disp('U')
U

end

   