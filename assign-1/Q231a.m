% 2.3.1 a part, 

disp('For the first case')
A = [3,1,1;1,7,3;2,0,4];
b = [5;11;5];
% The given matrix is invertible, moreover strictly diagonally dominant
% So both the Jacobi and the Gauss Seidel methods would converge
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

PJ = -inv(D)*(L + U);       % P Gauss Seidel
PGS = -inv(L+D)*U;                % P Jacobi 

disp('Spectral radius of P Gauss Seidel')
lambdamax_gs = max(abs(eig(PGS)))
disp('Spectral radius of P Jacobi')
lambdamaxj = max(abs(eig(PJ)))

% For the gauss Seidel P matrix, row(P) = 1, hence we can't say anything
% about convergence, whereas for the Jacobi method the spectral norm is less than 
% 1 hence it should converge for all intitial guesses. 

% For Jacobi and Gauss Seidel :: 
disp('The miminum number of iterations for 0.0001 relative error Gauss Seidel Method');
k = round(log(0.0001)/log(max(abs(eig(PGS)))))
disp('The miminum number of iterations for 0.0001 relative error Jacobi Method');
k = round(log(0.0001)/log(max(abs(eig(PJ)))))

fprintf('\n');

fprintf('\n');

fprintf('\n');
disp('Now for the second case')


A = [1,5,1;9,3,3;2,1,4];
b = [7;15;7];

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

PJ = -inv(D)*(L + U);       % P Gauss Seidel
PGS = -inv(L+D)*U;                % P Jacobi 

disp('Spectral radius of P Gauss Seidel')
lambdamax_gs = max(abs(eig(PGS)))
disp('Spectral radius of P Jacobi')
lambdamaxj = max(abs(eig(PJ)))


disp('For the Gauss Seidel and Jacobi Method we cannot infer anything as spectral radius of P is more than 1')

