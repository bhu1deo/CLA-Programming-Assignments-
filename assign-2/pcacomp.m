%% PCA Question CLA Assignment 

%% 1st part

I = im2double(imread('watch.bmp'));      % Read and convert to double 

R = I(:,:,1);
G = I(:,:,2);
B = I(:,:,3);

%% 2nd part, we will do for each channel separately 
[n1,n2] = size(R);
block_size = 8;            % 8x8 blocks 
num_blocks = n1*n2/(block_size^2);   % Number of Blocks for each channel

% R channel vectors 
numb_x = n1/block_size;
numb_y = n2/block_size;

i = 1;                       % Index in Matlab starts with 1 

R64 = zeros(block_size^2,1);       % The flattened vector for R channel 
count = 1;
while(i<=n1-block_size+1)
    j=1;
    while(j<=n2-block_size+1)
        R64(:,count) = reshape(R(i:i+block_size-1,j:j+block_size-1),1,[]);
        count = count + 1;  
        j = j + block_size; 
    end
    i = i + block_size;
end

% size(R64)             64 length vectors total 12288

% Repeat for B and G channels also
i = 1;
B64 = zeros(block_size^2,1);       % The flattened vector for B channel 
count = 1;
while(i<=n1-block_size+1)             % Extract blocks columnwise 
    j=1;
    while(j<=n2-block_size+1)
        B64(:,count) = reshape(B(i:i+block_size-1,j:j+block_size-1),1,[]);
        count = count + 1;  
        j = j + block_size; 
    end
    i = i + block_size;
end
% size(B64)

i = 1;
G64 = zeros(block_size^2,1);       % The flattened vector for G channel 
count = 1;
while(i<=n1-block_size+1)
    j=1;
    while(j<=n2-block_size+1)
        G64(:,count) = reshape(G(i:i+block_size-1,j:j+block_size-1),1,[]);
        count = count + 1;  
        j = j + block_size; 
    end
    i = i + block_size;
end
% size(G64)

%% 3rd part of the question 

sample_mean_R = sum(R64,2)/num_blocks;
sample_mean_G = sum(G64,2)/num_blocks;
sample_mean_B = sum(B64,2)/num_blocks;

%% 4th part of the Question 
% Do it for R-channel
Rcov = zeros(64,64);
for i=1:num_blocks
    Rcov = Rcov + (R64(:,i)-sample_mean_R)*transpose(R64(:,i)-sample_mean_R);
end
Rcov = Rcov/num_blocks; 
% Do it for G-channel
Gcov = zeros(64,64);
for i=1:num_blocks
    Gcov = Gcov + (G64(:,i)-sample_mean_G)*transpose(G64(:,i)-sample_mean_G);
end
Gcov = Gcov/num_blocks;
% Do it for B-channel
Bcov = zeros(64,64);
for i=1:num_blocks
    Bcov = Bcov + (B64(:,i)-sample_mean_B)*transpose(B64(:,i)-sample_mean_B);
end
Bcov = Bcov/num_blocks;

%% Sorting according to order

% R-channel 
[VR,LambdaR] = eig(Rcov);
LambdaR = diag(LambdaR); % Extract the diagonal elements
% Sort the eigenvalues in descending order, and get the indices.
[LambdaR,indices] = sort(LambdaR,'descend');
% Reorder the eigenvectors in the order given by the indices.
VR = VR(:,indices);

% G-channel 
[VG,LambdaG] = eig(Gcov);
LambdaG = diag(LambdaG); % Extract the diagonal elements
% Sort the eigenvalues in descending order, and get the indices.
[LambdaG,indices] = sort(LambdaG,'descend');
% Reorder the eigenvectors in the order given by the indices.
VG = VG(:,indices);


% B-channel 
[VB,LambdaB] = eig(Bcov);
LambdaB = diag(LambdaB); % Extract the diagonal elements
% Sort the eigenvalues in descending order, and get the indices.
[LambdaB,indices] = sort(LambdaB,'descend');
% Reorder the eigenvectors in the order given by the indices.
VB = VB(:,indices);

%% Image reconstruction using lesser number of principle components 

% Doing it for K = 5
% Doing This for R channel vectors 
% First find the coefficient matrix
disp('K=5')
K = 5;
R_channel_mean_matrix = repelem(sample_mean_R,1,num_blocks);  % Replicate
R_channel_mean_sub = R64 - R_channel_mean_matrix;    % Mean subtracted data 
coeff = zeros(64,num_blocks);                 % each col rep wrt the basis
coefficients = transpose(transpose(R_channel_mean_sub)*VR(:,1:K));   % Coeff. 
coeff(1:K,:) = coefficients;
restored_data_R = R_channel_mean_matrix + VR*coeff;


% Doing this for B channel vectors
B_channel_mean_matrix = repelem(sample_mean_B,1,num_blocks);  % Replicate
B_channel_mean_sub = B64 - B_channel_mean_matrix;    % Mean subtracted data 
coeff = zeros(64,num_blocks);         % each col rep wrt the basis
coefficients = transpose(transpose(B_channel_mean_sub)*VB(:,1:K));   % Coeff. 
coeff(1:K,:) = coefficients;
restored_data_B = B_channel_mean_matrix + VB*coeff;

% Doing this for G channel vectors 
G_channel_mean_matrix = repelem(sample_mean_G,1,num_blocks);  % Replicate
G_channel_mean_sub = G64 - G_channel_mean_matrix;    % Mean subtracted data 
coeff = zeros(64,num_blocks);                 % each col rep wrt the basis
coefficients = transpose(transpose(G_channel_mean_sub)*VG(:,1:K));   % Coeff. 
coeff(1:K,:) = coefficients;
restored_data_G = G_channel_mean_matrix + VG*coeff;


% Putting the blocks together 

restored_image_R = zeros(768,1024);
i = 1;

while(i<=n1-block_size+1)             % Extract blocks columnwise 
    j=1;
    while(j<=n2-block_size+1)
        restored_image_R(i:i+7,j:j+7) = reshape(restored_data_R(:,(i-1)*16+ceil(j/8)),[8,8]);
%         if(i<10 && i>1)
%             disp((i-1)*16+ceil(j/8))
%         end
        j = j + block_size; 
    end
    i = i + block_size; 
end

restored_image_B = zeros(768,1024);
i = 1;

while(i<=n1-block_size+1)             % Extract blocks columnwise 
    j=1;
    while(j<=n2-block_size+1)
        restored_image_B(i:i+7,j:j+7) = reshape(restored_data_B(:,(i-1)*16+ceil(j/8)),[8,8]);
        j = j + block_size; 
    end
    i = i + block_size; 
end

restored_image_G = zeros(768,1024);
i = 1;

while(i<=n1-block_size+1)             % Extract blocks columnwise 
    j=1;
    while(j<=n2-block_size+1)
        restored_image_G(i:i+7,j:j+7) = reshape(restored_data_G(:,(i-1)*16+ceil(j/8)),[8,8]);
        j = j + block_size; 
    end
    i = i + block_size; 
end

restored_image = cat(3, restored_image_R, restored_image_G, restored_image_B);
imshow(restored_image);
title('K=5');
%%
% K = 10
disp('K=10')
K = 10;
R_channel_mean_matrix = repelem(sample_mean_R,1,num_blocks);  % Replicate
R_channel_mean_sub = R64 - R_channel_mean_matrix;    % Mean subtracted data 
coeff = zeros(64,num_blocks);                 % each col rep wrt the basis
coefficients = transpose(transpose(R_channel_mean_sub)*VR(:,1:K));   % Coeff. 
coeff(1:K,:) = coefficients;
restored_data_R = R_channel_mean_matrix + VR*coeff;


% Doing this for B channel vectors
B_channel_mean_matrix = repelem(sample_mean_B,1,num_blocks);  % Replicate
B_channel_mean_sub = B64 - B_channel_mean_matrix;    % Mean subtracted data 
coeff = zeros(64,num_blocks);         % each col rep wrt the basis
coefficients = transpose(transpose(B_channel_mean_sub)*VB(:,1:K));   % Coeff. 
coeff(1:K,:) = coefficients;
restored_data_B = B_channel_mean_matrix + VB*coeff;

% Doing this for G channel vectors 
G_channel_mean_matrix = repelem(sample_mean_G,1,num_blocks);  % Replicate
G_channel_mean_sub = G64 - G_channel_mean_matrix;    % Mean subtracted data 
coeff = zeros(64,num_blocks);                 % each col rep wrt the basis
coefficients = transpose(transpose(G_channel_mean_sub)*VG(:,1:K));   % Coeff. 
coeff(1:K,:) = coefficients;
restored_data_G = G_channel_mean_matrix + VG*coeff;


% Putting the blocks together 

restored_image_R = zeros(768,1024);
i = 1;

while(i<=n1-block_size+1)             % Extract blocks columnwise 
    j=1;
    while(j<=n2-block_size+1)
        restored_image_R(i:i+7,j:j+7) = reshape(restored_data_R(:,(i-1)*16+ceil(j/8)),[8,8]);
%         if(i<10 && i>1)
%             disp((i-1)*16+ceil(j/8))
%         end
        j = j + block_size; 
    end
    i = i + block_size; 
end

restored_image_B = zeros(768,1024);
i = 1;

while(i<=n1-block_size+1)             % Extract blocks columnwise 
    j=1;
    while(j<=n2-block_size+1)
        restored_image_B(i:i+7,j:j+7) = reshape(restored_data_B(:,(i-1)*16+ceil(j/8)),[8,8]);
        j = j + block_size; 
    end
    i = i + block_size; 
end

restored_image_G = zeros(768,1024);
i = 1;

while(i<=n1-block_size+1)             % Extract blocks columnwise 
    j=1;
    while(j<=n2-block_size+1)
        restored_image_G(i:i+7,j:j+7) = reshape(restored_data_G(:,(i-1)*16+ceil(j/8)),[8,8]);
        j = j + block_size; 
    end
    i = i + block_size; 
end

restored_image = cat(3, restored_image_R, restored_image_G, restored_image_B);
imshow(restored_image);
title('K=10');
%%
% K = 20
disp('K=20')
K = 20;
R_channel_mean_matrix = repelem(sample_mean_R,1,num_blocks);  % Replicate
R_channel_mean_sub = R64 - R_channel_mean_matrix;    % Mean subtracted data 
coeff = zeros(64,num_blocks);                 % each col rep wrt the basis
coefficients = transpose(transpose(R_channel_mean_sub)*VR(:,1:K));   % Coeff. 
coeff(1:K,:) = coefficients;
restored_data_R = R_channel_mean_matrix + VR*coeff;


% Doing this for B channel vectors
B_channel_mean_matrix = repelem(sample_mean_B,1,num_blocks);  % Replicate
B_channel_mean_sub = B64 - B_channel_mean_matrix;    % Mean subtracted data 
coeff = zeros(64,num_blocks);         % each col rep wrt the basis
coefficients = transpose(transpose(B_channel_mean_sub)*VB(:,1:K));   % Coeff. 
coeff(1:K,:) = coefficients;
restored_data_B = B_channel_mean_matrix + VB*coeff;

% Doing this for G channel vectors 
G_channel_mean_matrix = repelem(sample_mean_G,1,num_blocks);  % Replicate
G_channel_mean_sub = G64 - G_channel_mean_matrix;    % Mean subtracted data 
coeff = zeros(64,num_blocks);                 % each col rep wrt the basis
coefficients = transpose(transpose(G_channel_mean_sub)*VG(:,1:K));   % Coeff. 
coeff(1:K,:) = coefficients;
restored_data_G = G_channel_mean_matrix + VG*coeff;


% Putting the blocks together 

restored_image_R = zeros(768,1024);
i = 1;

while(i<=n1-block_size+1)             % Extract blocks columnwise 
    j=1;
    while(j<=n2-block_size+1)
        restored_image_R(i:i+7,j:j+7) = reshape(restored_data_R(:,(i-1)*16+ceil(j/8)),[8,8]);
%         if(i<10 && i>1)
%             disp((i-1)*16+ceil(j/8))
%         end
        j = j + block_size; 
    end
    i = i + block_size; 
end

restored_image_B = zeros(768,1024);
i = 1;

while(i<=n1-block_size+1)             % Extract blocks columnwise 
    j=1;
    while(j<=n2-block_size+1)
        restored_image_B(i:i+7,j:j+7) = reshape(restored_data_B(:,(i-1)*16+ceil(j/8)),[8,8]);
        j = j + block_size; 
    end
    i = i + block_size; 
end

restored_image_G = zeros(768,1024);
i = 1;

while(i<=n1-block_size+1)             % Extract blocks columnwise 
    j=1;
    while(j<=n2-block_size+1)
        restored_image_G(i:i+7,j:j+7) = reshape(restored_data_G(:,(i-1)*16+ceil(j/8)),[8,8]);
        j = j + block_size; 
    end
    i = i + block_size; 
end

restored_image = cat(3, restored_image_R, restored_image_G, restored_image_B);
imshow(restored_image);
title('K=20');
%%
% K = 64
disp('K=64')
K = 64;
R_channel_mean_matrix = repelem(sample_mean_R,1,num_blocks);  % Replicate
R_channel_mean_sub = R64 - R_channel_mean_matrix;    % Mean subtracted data 
coeff = zeros(64,num_blocks);                 % each col rep wrt the basis
coefficients = transpose(transpose(R_channel_mean_sub)*VR(:,1:K));   % Coeff. 
coeff(1:K,:) = coefficients;
restored_data_R = R_channel_mean_matrix + VR*coeff;


% Doing this for B channel vectors
B_channel_mean_matrix = repelem(sample_mean_B,1,num_blocks);  % Replicate
B_channel_mean_sub = B64 - B_channel_mean_matrix;    % Mean subtracted data 
coeff = zeros(64,num_blocks);         % each col rep wrt the basis
coefficients = transpose(transpose(B_channel_mean_sub)*VB(:,1:K));   % Coeff. 
coeff(1:K,:) = coefficients;
restored_data_B = B_channel_mean_matrix + VB*coeff;

% Doing this for G channel vectors 
G_channel_mean_matrix = repelem(sample_mean_G,1,num_blocks);  % Replicate
G_channel_mean_sub = G64 - G_channel_mean_matrix;    % Mean subtracted data 
coeff = zeros(64,num_blocks);                 % each col rep wrt the basis
coefficients = transpose(transpose(G_channel_mean_sub)*VG(:,1:K));   % Coeff. 
coeff(1:K,:) = coefficients;
restored_data_G = G_channel_mean_matrix + VG*coeff;


% Putting the blocks together 

restored_image_R = zeros(768,1024);
i = 1;

while(i<=n1-block_size+1)             % Extract blocks columnwise 
    j=1;
    while(j<=n2-block_size+1)
        restored_image_R(i:i+7,j:j+7) = reshape(restored_data_R(:,(i-1)*16+ceil(j/8)),[8,8]);
%         if(i<10 && i>1)
%             disp((i-1)*16+ceil(j/8))
%         end
        j = j + block_size; 
    end
    i = i + block_size; 
end

restored_image_B = zeros(768,1024);
i = 1;

while(i<=n1-block_size+1)             % Extract blocks columnwise 
    j=1;
    while(j<=n2-block_size+1)
        restored_image_B(i:i+7,j:j+7) = reshape(restored_data_B(:,(i-1)*16+ceil(j/8)),[8,8]);
        j = j + block_size; 
    end
    i = i + block_size; 
end

restored_image_G = zeros(768,1024);
i = 1;

while(i<=n1-block_size+1)             % Extract blocks columnwise 
    j=1;
    while(j<=n2-block_size+1)
        restored_image_G(i:i+7,j:j+7) = reshape(restored_data_G(:,(i-1)*16+ceil(j/8)),[8,8]);
        j = j + block_size; 
    end
    i = i + block_size; 
end

restored_image = cat(3, restored_image_R, restored_image_G, restored_image_B);
imshow(restored_image);
title('K=64');
%% Plotting the error trajectory for different values of K








