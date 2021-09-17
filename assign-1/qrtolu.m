% QR to LU without explicitly constructing the original matrix
% Find QR and then find LU of Q, you are done
% product of upper triangular matrices is a upper triangular matrix
% Find the inverse and the determinant 

% We do the following::
% 1.) First find the QR decomposition 
% 2.) Find LU decomposition of Q, i.e PQ = LU
% 3.) Thus the final form becomes A = QR = PTLUR = PT(LU')
% 4.) Which is equivalent to PA = LU', here U' is also an upper triangular
% matrix and the problem is complete 

% Matrix A


function [L,U,P] = qrtolu(Q,R)        % Pass QR get L,U,P; don't touch A
format 

% LU of Q

[determinant,inverse,L,U,P] = LUpartial(Q);

U = U*R;            % L and P remain same 

disp('Q*R')
Q*R
disp('PTLUR')
transpose(P)*L*U*R
disp('Error NORM')           % See the norm of the error 
norm(transpose(P)*L*U-Q*R)

end