% LU with partial pivoting 
% Return the determinant and inverse of the matrix 
% The input matrix is square and invertible 
% In each iteration, choose the pivot to be the largest value in that
% column interchange rows if necessary,keep track of the permutations 
% The LU of PA needs to be carried out where P is the
% permutation(transposition) matrix 

function [determinant,inverse,L,U,P] = LUpartial(A)     % input is an invertible matrix
format 

if(det(A)==0)           % Matrix is nonsingular 
    disp('LU not possible, matrix not invertible')
    return;
end
% If matrix is non-singular then LU with pivoting is possible
interchange = 0;
n1 = size(A);
n = n1(1);                               % Rows = Columns here
C = A;
n1 = size(A);
n = n1(1);                               % Rows = Columns here
% first find the first pivot in the first column 
[val,idx] = max(abs(A(:,1))); 

if(idx~=1)
    interchange = interchange + 1;
end
% interchange the rows 
tmp = A(1,:);
A(1,:) = A(idx,:);
A(idx,:) = tmp;
P = eye(n);                 % The permutation matrix
tmp = P(1,:);
P(1,:) = P(idx,:);
P(idx,:) = tmp;        
% try to work with multipliers here instead
L = zeros(n); 
U = zeros(n);
L(2:n,1) = A(2:n,1)/A(1,1);
U(1,:) = A(1,:); 
for j = 2:n
   A(j,:) = A(j,:)-L(j,1)*A(1,:);
   A(j,1) = 0;
end
B = A(2:n,2:n);          % Operate on the submatrix now     
for k=2:n
    s = size(B);
    if(s(1)==1 && s(2)==1)            % Last submatrix
        U(n,n) = B;            % shouldn't be 0 
        break
    else
        % First Do an interchange 
        [val,idx] = max(abs(B(:,1)));  
        if(idx~=1)
            interchange = interchange + 1;
        end
        % interchange the rows 
        tmp = B(1,:);
        B(1,:) = B(idx,:);
        B(idx,:) = tmp;
        pidx = idx+k-1;
        permut = eye(n);
        tmp = permut(k,:);         % interchanging wrt the kth row 
        permut(k,:) = permut(pidx,:);
        permut(pidx,:) = tmp;
        P = permut*P;
        tmp = L(k,:);
        L(k,:) = L(pidx,:);
        L(pidx,:) = tmp;
        L(k+1:n,k) = B(2:s(1),1)/B(1,1);
        
        U(k,k:n) = B(1,:);
        for j = 2:s(1)
            B(j,:) = B(j,:)-L(j+k-1,k)*B(1,:);
            B(j,1) = 0;
        end
        
    end
    
    B = B(2:s(1),2:s(1));
end
L = L + eye(n);
A = C;
% Following things are for display purposes 
% disp("P*A");
% P*A 
% disp("L*U");
% L*U
% disp("A");
% A
% disp("PTLU");
% transpose(P)*L*U             % The original matrix 
% disp("NORM PA-LU");
% norm(P*A-L*U)

% Now compute the determinant of A 
if(mod(interchange,2)==1)
    determinantP = -1;
else
    determinantP = 1;
end
determinant = determinantP*1*prod(diag(U));

disp(determinant);

% Now solve the inverse problem 
% LUx = ei ---> Ly = ei & Ux = y  

% For each ei solve for y:: 

inverse_mat = zeros(n); 
T = eye(n);
for j = 1:n
    t = T(:,j);              % this is ej
    for i = 1:n
       if(i==1)
           y = t(1,:);
       else
           y(i,:) = t(i,:) - (L(i,1:i-1))*(y);     % Till i-1 we have computed 
       end
    end
    % Now we have yi, we find xi and then compose the matrix
    
    for k = 1:n
       id = n+1-k; 
       if(k==1)
           x = y(id,:)/U(id,id);
       else
           xtemp = (y(id,:) - (U(id,id+1:n))*(x))/U(id,id);
           x = vertcat(xtemp,x);
       end
    end
    
 

    inverse_mat(:,j) = x;      
end
inverse = inverse_mat*P;

end






