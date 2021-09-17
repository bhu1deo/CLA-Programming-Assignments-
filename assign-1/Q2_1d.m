% Q2.1d part

n = 3;

% Generate 2 orthonormal vectors

u1 = rand(n,1);
u2 = rand(n,1);

v1 = u1/norm(u1);
r2 = u2 - dot(u2,v1)*v1;
v2 = r2/norm(r2);

% Verify Orthonormality of v1,v2
disp('The norm of v1 is')
norm(v1)
disp('The norm of v2 is')
norm(v2)
disp('The inner product between v1 and v2 is')
dot(v1,v2)

% Matrix A
A = 50000*v1*transpose(v1)+2*v2*transpose(v2);

% x and b

b = A*rand(n,1);


[Qgs,Rgs] = gs(A);
[Qhr,Rhr] = hr(A);
[Qprop,Rprop] = prop(A);
% Solve Rx = QTb
bgs = transpose(Qgs)*b;
bhr = transpose(Qhr)*b;
bprop = transpose(Qprop)*b;

% Estimate of x by Gram Schmidt
for k=1:n
    id = n+1-k;           % n to 1 reverse back-substitution 
    if(k==1)
        xgs = bgs(id,:)/Rgs(id,id);
    else
        xtemp = (bgs(id,:)-(Rgs(id,id+1:end))*(xgs))/Rgs(id,id);
        xgs = vertcat(xtemp,xgs);
    end
end


% Estimate of x by Reflection
for k=1:n
    id = n+1-k;           % n to 1 reverse back-substitution 
    if(k==1)
        xhr = bhr(id,:)/Rhr(id,id);
    else
        xtemp = (bhr(id,:)-(Rhr(id,id+1:end))*(xhr))/Rhr(id,id);
        xhr = vertcat(xtemp,xhr);
    end
end


% Estimate of x by Rotation (PROP)
for k=1:n
    id = n+1-k;           % n to 1 reverse back-substitution 
    if(k==1)
        xprop = bprop(id,:)/Rprop(id,id);
    else
        xtemp = (bprop(id,:)-(Rprop(id,id+1:end))*(xprop))/Rprop(id,id);
        xprop = vertcat(xtemp,xprop);
    end
end
% Rotations method is the best to solve Ax=b, followed by slightly worse 
% performance of reflections and BAD performance by Gram Schmidt
disp('Ax-b for Gram Schmidt')
norm(A*xgs-b)
disp('Ax-b for Reflections')
norm(A*xhr-b)
disp('Ax-b for Rotations')
norm(A*xprop-b)

