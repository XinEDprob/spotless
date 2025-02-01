%% Test for minimum eigenvalue of M-tensor

%
% cp = [2, -1, 0, -2, 3];
% cp = [2.5, -2.5, 0, -0.5, 1.5];
% cp = [5, -2.5/3, 0, -1/3, 4.5];

% exemple 1

% symmetric
A(:,:,1,1) = [3, -1/2;
           -1/2, 0];
A(:,:,1,2) = [-1/2, 0;
           0, -1/4];
A(:,:,2,1) = [-1/2, 0;
           0, -1/4];
A(:,:,2,2) = [0, -1/4;
           -1/4, 2];
% 
B(:,:,1,1) = [1.5, -0.5/4;
           -0.5/4, 0];
B(:,:,1,2) = [-0.5/4, 0;
           0, -2.5/4];
B(:,:,2,1) = [-0.5/4, 0;
           0, -2.5/4];
B(:,:,2,2) = [0, -2.5/4;
           -2.5/4, 2.5];
% 
tauA = 0.97225;
tauB = 0.5;
n = 2;       
       
% asymmetric       
% A(:,:,1,1) = [3, 0;
%            -2, 0];
% A(:,:,1,2) = [0, 0;
%            0, 0];
% A(:,:,2,1) = [0, 0;
%            0, 0];
% A(:,:,2,2) = [0, -1;
%            0, 2];
% 
% B(:,:,1,1) = [1.5, 0;
%            -0.5, 0];
% B(:,:,1,2) = [0, 0;
%            0, 0];
% B(:,:,2,1) = [0, 0;
%            0, 0];
% B(:,:,2,2) = [0, -1;
%            -1.5, 2.5];
%        
%        
% tauA = 1;
% tauB = 0.5;
% n = 2;

%%% example 2
% symmetric
% A(:,:,1,1) = [3.8, -0.5;
%            -0.5, -13/30];
% A(:,:,1,2) = [-0.5, -13/30;
%            -13/30, -0.5];
% A(:,:,2,1) = [-0.5, -13/30;
%            -13/30, -0.5];
% A(:,:,2,2) = [-13/30, -0.5;
%            -0.5, 3.9];
% 
% B(:,:,1,1) = [3.2, -0.675;
%            -0.675, -1/3];
% B(:,:,1,2) = [-0.675, -1/3;
%            -1/3, -0.35];
% B(:,:,2,1) = [-0.675, -1/3;
%            -1/3, -0.35];
% B(:,:,2,2) = [-1/3, -0.35;
%            -0.35, 3.9];
% tauA = 0.54995;
% tauB = 0.41253;
% n = 2;

% asymmetric
% A(:,:,1,1) = [3.8, -0.8;
%            -0.4, -0.8];
% A(:,:,1,2) = [-0.7, -0.3;
%            -0.6, -0.4];
% A(:,:,2,1) = [-0.1, -0.5;
%            -0.4, -0.4];
% A(:,:,2,2) = [0, -0.2;
%            -1, 3.9];
% 
% B(:,:,1,1) = [3.2, -0.2;
%            -0.6, -0.2];
% B(:,:,1,2) = [-1, -0.5;
%            -0.7, -0.5];
% B(:,:,2,1) = [-0.9, -0.5;
%            0, -0.2];
% B(:,:,2,2) = [-0.1, 0;
%            -0.7, 3.9];
% %    
% % 
% tauA = 0.6525;
% tauB = 0.4158;
% n = 2;


%%example 3

% A(:,:,1,1) = [3.8, -0.575;
%            -0.575, -11/30];
% A(:,:,1,2) = [-0.575, -11/30;
%            -11/30, -0.4];
% A(:,:,2,1) = [-0.575, -11/30;
%            -11/30, -0.4];
% A(:,:,2,2) = [-11/30, -0.4;
%            -0.4, 3.7];
% 
% B(:,:,1,1) = [3.5, -0.35;
%            -0.35, -23/60];
% B(:,:,1,2) = [-0.35, -23/60;
%            -23/60, -0.525];
% B(:,:,2,1) = [-0.35, -23/60;
%            -23/60, -0.525];
% B(:,:,2,2) = [-23/60, -0.525;
%            -0.525, 3.1];
% tauA = 0.69696;
% tauB = 0.3717;
% n = 2;

% asymmetric
% A(:,:,1,1) = [3.8, -0.6;
%               -0.9, -0.4];
% A(:,:,1,2) = [-0.6, -0.4;
%            -0.2, -0.6];
% A(:,:,2,1) = [-0.2, -0.5;
%            -0.4, -0.1];
% A(:,:,2,2) = [-0.3, -0.6;
%            -0.3, 3.7];
% 
% B(:,:,1,1) = [3.5, -0.1;
%            -0.3, -0.6];
% B(:,:,1,2) = [-0.2, -0.4;
%            -0.3, -0.5];
% B(:,:,2,1) = [-0.8, -0.1;
%            -0.4, -0.3];
% B(:,:,2,2) = [-0.5, -0.5;
%            -0.9, 3.1];
%        
% tauA = 0.6972;
% tauB = 0.3234;
% n = 2;
%% bounds 3.1
% lower bound

min_1 = Inf;
for i = 1:n
    min_1 = min(min_1, A(i,i,i,i)*tauB + B(i,i,i,i)*tauA);
%     fprintf('min: %d', min_1)
end
min_1 = min_1 - tauA*tauB;
fprintf('bounds 3.1\n')
fprintf('min: %d\n',min_1);

%% 


max_sliceA = slice_MAX(A, n);

max_sliceB = slice_MAX(B, n);

min_sliceA = slice_MIN(A, n);

min_sliceB = slice_MIN(B, n);

gamma1_AB = GAMMA1(A,B,max_sliceA,tauB,n);

gamma1_BA = GAMMA1(B,A,max_sliceB,tauA,n);

gamma2_AB = GAMMA2(A,B,min_sliceA,tauB,n);

gamma2_BA = GAMMA2(B,A,min_sliceB,tauA,n);

gamma3 = GAMMA3(A,B,max_sliceA, max_sliceB, tauA, tauB,n);

gamma4 = GAMMA4(A,B,min_sliceA, min_sliceB, tauA, tauB,n);


%% bounds 3.3
fprintf('bounds 3.3\n')
fprintf('min: %d\n',gamma1_AB);
fprintf('max: %d\n',gamma2_AB);

%% bounds 3.4
fprintf('bounds 3.4\n')
fprintf('min: %d\n',gamma1_BA);
fprintf('max: %d\n',gamma2_BA);

%% bounds 3.5
fprintf('bounds 3.5\n')
fprintf('min: %d\n',max(gamma1_AB, gamma1_BA));
fprintf('max: %d\n',min(gamma2_AB, gamma2_BA));

%% bounds 3.6
fprintf('bounds 3.6\n')
fprintf('min: %d\n',gamma3);
fprintf('max: %d\n',gamma4);

%% bounds 3.7
min_1 = Inf;
slice_sum_A = slice_SUM(A, n);
slice_sum_B = slice_SUM(B, n);
for i = 1:n
    min_1 = min(min_1, A(i,i,i,i)*B(i,i,i,i) - slice_sum_A(i)*slice_sum_B(i));
end
fprintf('bounds 3.7 \n')
fprintf('min: %d\n',min_1);

function max_slices =  slice_MAX(A, n)
    max_slices = zeros(n,1);
    for i = 1:n
        max_slices(i) = -Inf;
        for i2 = 1:n
            for i3 = 1:n
                for i4 = 1:n
                    if ~(i2== i && i3==i && i4== i)
                        max_slices(i) = max(max_slices(i), abs(A(i,i2,i3,i4)));
                    end
                end
            end
        end
    end
end

function min_slices =  slice_MIN(A, n)
    min_slices = zeros(n,1);
    for i = 1:n
        min_slices(i) = Inf;
        for i2 = 1:n
            for i3 = 1:n
                for i4 = 1:n
                    if ~(i2== i && i3==i && i4== i)
                        min_slices(i) = min(min_slices(i), abs(A(i,i2,i3,i4)));
                    end
                end
            end
        end
    end
end

function gamma1 = GAMMA1(A, B, max_sliceA, tauB, n)
    gamma1 = Inf;
%     max_sliceA = slice_MAX(A);
    for i = 1:n
        gamma1 = min(gamma1, (A(i,i,i,i) - max_sliceA(i))*B(i,i,i,i) + max_sliceA(i)*tauB);
    end
end

function gamma2 = GAMMA2(A, B, min_sliceA, tauB, n)
    gamma2 = -Inf;
%     max_sliceA = slice_MAX(A);
    for i = 1:n
        gamma2 = max(gamma2, (A(i,i,i,i) - min_sliceA(i))*B(i,i,i,i) + min_sliceA(i)*tauB);
    end
end


function gamma3 = GAMMA3(A, B, max_sliceA, max_sliceB, tauA, tauB, n)
    gamma3 = Inf;
%     max_sliceA = slice_MAX(A);
    for i = 1:n
        gamma3 = min(gamma3, (A(i,i,i,i)*B(i,i,i,i) - sqrt(max_sliceA(i))*sqrt(max_sliceB(i))*...
            sqrt(A(i,i,i,i) - tauA)*sqrt(B(i,i,i,i) - tauB)));
    end
end


function gamma4 = GAMMA4(A, B, min_sliceA, min_sliceB, tauA, tauB, n)
    gamma4 = -Inf;
%     max_sliceA = slice_MAX(A);
    for i = 1:n
        gamma4 = max(gamma4, (A(i,i,i,i)*B(i,i,i,i) - sqrt(min_sliceA(i))*sqrt(min_sliceB(i))*...
            sqrt(A(i,i,i,i) - tauA)*sqrt(B(i,i,i,i) - tauB)));
    end
end

function sum_slices =  slice_SUM(A, n)
    sum_slices = zeros(n,1);
    for i = 1:n
        for i2 = 1:n
            for i3 = 1:n
                for i4 = 1:n
                    if ~(i2== i && i3==i && i4== i)
                        sum_slices(i) = sum_slices(i) +  A(i,i2,i3,i4);
                    end
                end
            end
        end
    end
end