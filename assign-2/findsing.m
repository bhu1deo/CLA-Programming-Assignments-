% Now SVD of a general rectangular matrix
% The steps are as are as follows:
% 1.) Convert to bidiagonal form
% 2.) Use QR Algorithm (basic or with shifts) to get to singular Values 
% We will follow the Golub-Reinsch method and avoid complications 
% Singular values of bidiagonal matrix are same as that of A

% Also we use rotation matrices again for QR decomposition as bidiag matrix
% is sparse.
%% We will do here for tall matrices, FAT matrices can be easily generalized 

function [singular_values] = findsing(A)
    % First convert to bidiagonal form by multiplying left and right with
    % orthogonal singular values is invariant to such manipulations here
    % assume m>n here
    
    siz = size(A);
    B = A;
    m = siz(1);           % Number of rows
    n = siz(2);           % Number of columns m>n
    % We will NOT store the U,V vectors as just singular values are needed
    k = 1;     % The row index A is an mxn matrix
    
    while(k<n-1)       % Due to bidiagonality
        % kill the column of A 
        x = A(k:end,k);   % kth column of A
        xsiz = size(x);
        for j = 2:xsiz(1)
            Q = eye(m);
            Q(k,k) = x(1)/sqrt(x(1)^2+x(j)^2);
            Q(k-1+j,k-1+j) = x(1)/sqrt(x(1)^2+x(j)^2);
            Q(k-1+j,k) = -x(j)/sqrt(x(1)^2+x(j)^2);
            Q(k,k-1+j) = x(j)/sqrt(x(1)^2+x(j)^2);
            A = Q*A;
            x = A(k:end,k); 
        end
        
        % kill the row of A
        x = A(k,k+1:end);
        xsiz = size(x);
        if(xsiz(2)~=2)
            for j = 2:xsiz(2)
                Q = eye(n);
                Q(k+1,k+1) = x(1)/sqrt(x(1)^2+x(j)^2);
                Q(k+j,k+j) = x(1)/sqrt(x(1)^2+x(j)^2);
                Q(k+1,k+j) = -x(j)/sqrt(x(1)^2+x(j)^2);
                Q(k+j,k+1) = x(j)/sqrt(x(1)^2+x(j)^2);
                A = A*Q;
                x = A(k,k+1:end);
            end
        else
           % Last row kill it and also kill the remaining columns 
           for j = 2:xsiz(2)
                Q = eye(n);
                Q(k+1,k+1) = x(1)/sqrt(x(1)^2+x(j)^2);
                Q(k+j,k+j) = x(1)/sqrt(x(1)^2+x(j)^2);
                Q(k+1,k+j) = -x(j)/sqrt(x(1)^2+x(j)^2);
                Q(k+j,k+1) = x(j)/sqrt(x(1)^2+x(j)^2);
                A = A*Q
                x = A(k,k+1:end);
           end
           % kill remaining columns
     
           for mi = k+1:n
               x = A(mi:end,mi); 
               xsiz = size(x);
               for j = 2:xsiz(1)
                    Q = eye(m);
                    Q(mi,mi) = x(1)/sqrt(x(1)^2+x(j)^2);
                    Q(mi-1+j,mi-1+j) = x(1)/sqrt(x(1)^2+x(j)^2);
                    Q(mi,mi-1+j) = x(j)/sqrt(x(1)^2+x(j)^2);
                    Q(mi-1+j,mi) = -x(j)/sqrt(x(1)^2+x(j)^2);
                    A = Q*A;
                    x = A(mi:end,mi); 
               end 
           end
           break
           
           
        end
        
        k = k + 1;    
    end
    bidiag = A;
    disp('The bidiagonal equivalent is');
%     bidiag
   
    tri_diag = transpose(bidiag)*bidiag;
    
    % Now Apply QR Algorithm on tridiagonal form to get singular values 
    for h=1:1000    % This needs to be increased for large n 
        [Q,R] = QR_tridiag(tri_diag);
        tri_diag = R*Q;
    end
    disp('The Singular Values of A are')
    singular_values = sqrt(diag(tri_diag))
end
