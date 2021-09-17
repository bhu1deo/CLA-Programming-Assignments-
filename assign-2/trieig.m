% Write a program (trieig.m) to implement QR algorithm and thereby 
% find eigenvalues of tridiagonal matrix. 
% Note that QR decomposition should be performed using Householder trans-
% formation or method proposed in PGA1. 
% Also the subroutine function to implement QR
% decomposition should be tailor-made for tridiagonal matrix. 
% Print per iteration complexity for QR algorithm in this specific case.


%% We will implement the basic QR algorithm to find the eigenvalues 
% For this specific case of tridiagonal matrices we use given's rotators 
% We will require n-1 ortho matrices for the complete matrix
function [eigenvalues] = trieig(A)
% The input matrix is assumed to be symmetric tri-diagonal matrix 
% A = QR then A = RQ 
disp('Note : Ensure Order of found eigenvalues is correct before computing error norm');
fprintf('\n');
for h=1:1000
   [Q,R] = QR_tridiag(A);
   A = R*Q;
end
eigenvalues = diag(A);  % Directly the diagonal elements
disp('The per iteration complexity for A=QR is O(n^2)')
fprintf('\n');
disp('This is because we use only n-1 total givens rotators for every iteration')
disp('Computing R*Q is then also O(n^2) as R is upper triangular')


end