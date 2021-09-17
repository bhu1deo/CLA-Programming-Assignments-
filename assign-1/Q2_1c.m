% Q2.1c part

n = 40;
A = 0.00001*eye(n)+hilb(n);

% The 3 methods 
[Qgs,Rgs] = gs(A);
[Qhr,Rhr] = hr(A);
[Qprop,Rprop] = prop(A);

% Difference between QTQ and Identity

disp('Error in QTQ and I for Gram Schmidt');      % Large Error 
norm(transpose(Qgs)*Qgs-eye(n))
disp('Error in QTQ and I for Reflection');           % Vastly Better
norm(transpose(Qhr)*Qhr-eye(n))
disp('Error in QTQ and I for Rotation');             % The Best
norm(transpose(Qprop)*Qprop-eye(n))


% Error between A and QR::
disp('Error between A and QR Gram Schmidt');
norm(A-Qgs*Rgs)
disp('Error between A and QR Reflection');
norm(A-Qhr*Rhr)
disp('Error between A and QR Rotation');
norm(A-Qprop*Rprop)

