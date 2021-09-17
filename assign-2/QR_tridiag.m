% QR for tridiagonal matrices 
% Use n-1 givens rotators for nxn matrix specifically for tridiag matrices


function [Q,R] = QR_tridiag(A)
% A is a symmetric tridiagonal matrix 
n = size(A);
n = n(1);                       % Rows = Columns here 
Q = eye(n);
R = A;    % working matrix
for k=1:n-1
    % Do just k-k+1 plane rotation Also rotation preserves NORM
    Qstar = eye(n);                           
    x = R(k:end,k);
    c = x(1,:)/sqrt(x(1,:)^2+x(2,:)^2);
    s = x(2,:)/sqrt(x(1,:)^2+x(2,:)^2);
    Qstar(k,k) = c;
    Qstar(k+1,k+1) = c;
    Qstar(k,k+1) = s;
    Qstar(k+1,k) = -s;
    Q = Qstar*Q;
    R = Qstar*R;
end
Q = transpose(Q);

end
