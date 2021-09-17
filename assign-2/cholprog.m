% Cholesky Decomposition 

function [L] = cholprog(A)
% A is PD symmetric square and real; it's assumed not checked
n = size(A);
n = n(1);
L = zeros(n);

for k=1:n
   % compute the diagonal entry first 
   if(k==1)
      L(1,1) = sqrt(A(1,1));        % Easy here 
   else
      
      L(k,k) = sqrt(A(k,k)-norm(L(k,1:k-1))^2);   % Similar to backsubs
   end
   % Now L(k,k) is found
   % Now find the columns of L
   
   if(k==n)                 % Only L(n,n) is enough here 
       break
   end
   
   if(k==1)
      L(2:end,1) = A(1,2:end)/L(1,1);   % First col. of L is found
   else
       for i = k+1:n
%            
           L(i,k) = (A(i,k)-dot(L(i,1:k-1),L(k,1:k-1)))/L(k,k);
       end
       
          
       % All previous cols of L are found out; we have to use this fact
   end
   
   
    
end

end