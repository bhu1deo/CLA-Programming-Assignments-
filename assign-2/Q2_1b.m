% CLA Programming Assignment-2 Q2_1b.m 

%% We will define a local function here to aid the process The error traj is plotted

A = [3,3,4;3,7,6;4,6,10];
larg_eig = comp_largest_eigenvalue(A);
function [largest_eig] = comp_largest_eigenvalue(A)
    % We will plot the largest eigenvalue norm error here
    max_eig = max(eig(A)); % computed by Matlab
    error_norm = 0;
    % convert it to tridiagonal form
    Z = sym2tri(A);  % Tri-Diagonal Form of A
    
    for h=1:10
        [Q,R] = QR_tridiag(Z);
        Z = R*Q;
        eig_max = Z(1,1);
        error_norm(end+1) = norm(eig_max-max_eig);
    end
    largest_eig = eig_max;
    error_norm = error_norm(2:end);
    X = 1:1:10;
    plot(X,error_norm)
end