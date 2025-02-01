
%% Table 1-3 (our methods) 

N = 3;
m = 3;

% The followings are two examples, to test one of them, you just need to
% comment out another one
% symmetrized Example 3.1 from https://www.sciencedirect.com/science/article/pii/S002437951830404X
% T1(1,:,:) = [7, -1/3, -1/3;
%     -1/3, -2.1, -7/6;
%     -1/3, -7/6, -2/3];
% T1(2,:,:) = [-1/3, -2.1, -7/6;
%     -2.1, 12, -1/3;
%     -7/6, -1/3, -6.5/3];
% T1(3,:,:) = [-1/3, -7/6, -2/3;
%     -7/6, -1/3, -6.5/3;
%     -2/3, -6.5/3, 50];

% symmetrized Example 3.2 from https://www.sciencedirect.com/science/article/pii/S002437951830404X
T1(1,:,:) = [12, -0.5/3, -0.1;
     -0.5/3, -1.6, -2.5;
    -0.1, -2.5, -3.5/3];
T1(2,:,:) = [ -0.5/3, -1.6, -2.5;
    -1.6, 30, -1/3;
    -2.5, -1/3, -7/3];
T1(3,:,:) = [-0.1, -2.5, -3.5/3;
    -2.5, -1, -7/3;
    -3.5/3, -7/3, 15];



n = N;
n_var = nchoosek(n+m-1, m);
Qt = symtensor([1: n_var], m, n);

unique_monomial = [];

Qt_inv = containers.Map;
Qt_map = containers.Map;
i = 1;
for i1 = 1:n
    for i2 = i1:n
        for i3 = i2:n
            Qt_inv(int2str(i)) = [i1,i2,i3];
%             Qt_map()
            unique_monomial = [unique_monomial; [i1,i2,i3]];
            i = i+1;
        end
    end
end

polycoeff = zeros(n_var, 1);
for i1 = 1:n
    for i2 = 1:n
        for i3 = 1:n
            for j = 1:size(unique_monomial, 1)
                if sort([i1, i2, i3]) == unique_monomial(j,:)
                    polycoeff(j) = polycoeff(j) + T1(i1,i2,i3);
                end
            end
        end
    end
end


% [sol, time] = SDDT(cp, N, m);

[sol, time] = fast_SDDT33(polycoeff);

%% Example 4
cp = zeros(1, 20);
cp(6) = -3;
cp(9) = -3;
m = 3;
N = 4;

[sol, time] = SDDT(cp, N, m);

%% Example 5
m = 3;
n = 2;

% empty tensors
D1 = zeros(n,n,n);
for i1 = 1:n
    for i2 = 1:n
        for i3 = 1:n
            D1(i1, i2, i3) = abs(tan(i1+i2+i3)) ;
        end
    end
end

% diagonal term
s1 = 0;
for i = 1:n
    s1 = max(s1, sum(D1(i,:,:),'all'));
end

s1 = s1*1.01;
% s1 = 1500;

% Matrix
T1 = zeros(n,n,n);
for i1 = 1:n
    for i2 = 1:n
        for i3 = 1:n
            if i1 == i2 && i2== i3
                T1(i1, i2, i3) = s1 - abs(tan(i1+i2+i3));
            else
               T1(i1, i2, i3) = -abs(tan(i1+i2+i3)) ; 
            end
        end
    end
end


n_var = nchoosek(n+m-1, m);
Qt = symtensor([1: n_var], m, n);

unique_monomial = [];

Qt_inv = containers.Map;
Qt_map = containers.Map;
i = 1;
for i1 = 1:n
    for i2 = i1:n
        for i3 = i2:n
            Qt_inv(int2str(i)) = [i1,i2,i3];
%             Qt_map()
            unique_monomial = [unique_monomial; [i1,i2,i3]];
            i = i+1;
        end
    end
end

polycoeff = zeros(n_var, 1);
for i1 = 1:n
    for i2 = 1:n
        for i3 = 1:n
            for j = 1:size(unique_monomial, 1)
                if sort([i1, i2, i3]) == unique_monomial(j,:)
                    polycoeff(j) = polycoeff(j) + T1(i1,i2,i3);
                end
            end
        end
    end
end

[sol, time] = SDDT(polycoeff, n, m);

%% Example 6
m = 3;
n = 3;

h = 1/(n-1);

% empty tensors
T1 = zeros(n,n,n);

T1(1,1,1) = 1/(h*h);
T1(n,n,n) = 1/(h*h);

for i = 2:n-1
    T1(i,i,i) = 2/(h*h);
end

for i = 2:n-1
    T1(i,i-1,i) = -1/(h*h*(m-1));
    T1(i,i,i-1) = -1/(h*h*(m-1));
    T1(i,i+1,i) = -1/(h*h*(m-1));
    T1(i,i,i+1) = -1/(h*h*(m-1));
end

n_var = nchoosek(n+m-1, m);
Qt = symtensor([1: n_var], m, n);

polycoeff = [4, 0, 0, -4, 0, 0, 8, -4, 0 ,4];

% a symmetrillization tensor of T1
T11 = zeros(n,n,n);
T11(1,1,1) = 1/(h*h);
T11(n,n,n) = 1/(h*h);

for i = 2:n-1
    T11(i,i,i) = 2/(h*h);
end
T11(1,2,2) = -4/3;
T11(2,1,2) = -4/3;
T11(2,2,1) = -4/3;
T11(2,2,3) = -4/3;
T11(2,3,2) = -4/3;
T11(3,2,2) = -4/3;

[sol, time] = SDDT(polycoeff, n, m);

%% Example 7
m = 3;
n = 3;

beta = 0.3;

% The followings are two examples, to test one of them, you just need to
% comment out another one
% the first example of Example 7
% matrix1 = [0.5810, 0.2432, 0.1429;
%            0, 0.4109, 0.0701;
%            0.4190, 0.3459, 0.7870];
% 
% matrix2 = [0.4708, 0.1330, 0.0327;
%            0.1341, 0.5450, 0.2042;
%            0.3951, 0.3220, 0.7631];
% 
% matrix3 = [0.4381, 0.1003, 0;
%            0.0229, 0.4338, 0.0930;
%            0.5390, 0.4659, 0.9070;];

% the second example in Example 7
matrix1 = [0.9000, 0.3340, 0.3106;
           0.0690, 0.6108, 0.0754;
           0.0310, 0.0552, 0.6140];

matrix2 = [0.6700, 0.1040, 0.0805;
           0.2892, 0.8310, 0.2956;
           0.0408, 0.0650, 0.6239];

matrix3 = [0.6604, 0.0945, 0.0710
           0.0716, 0.6133, 0.0780;
           0.2680, 0.2922, 0.8501;];

% Create the 3D tensor P
P = cat(3, matrix1, matrix2, matrix3);

% empty tensors
D1 = zeros(n,n,n);
for i1 = 1:n
    for i2 = 1:n
        for i3 = 1:n
            D1(i1, i2, i3) = beta*P(i1, i2, i3) ;
        end
    end
end

% diagonal term
s1 = 0;
for i = 1:n
    s1 = max(s1, sum(D1(i,:,:),'all'));
end

s1 = 1;

% Matrix
T1 = zeros(n,n,n);
for i1 = 1:n
    for i2 = 1:n
        for i3 = 1:n
            if i1 == i2 && i2== i3
                T1(i1, i2, i3) = s1 - D1(i1, i2, i3);
            else
               T1(i1, i2, i3) = - D1(i1, i2, i3); 
            end
        end
    end
end


n_var = nchoosek(n+m-1, m);
Qt = symtensor([1: n_var], m, n);

unique_monomial = [];

Qt_inv = containers.Map;
Qt_map = containers.Map;
i = 1;
for i1 = 1:n
    for i2 = i1:n
        for i3 = i2:n
            Qt_inv(int2str(i)) = [i1,i2,i3];
%             Qt_map()
            unique_monomial = [unique_monomial; [i1,i2,i3]];
            i = i+1;
        end
    end
end

polycoeff = zeros(n_var, 1);
for i1 = 1:n
    for i2 = 1:n
        for i3 = 1:n
            for j = 1:size(unique_monomial, 1)
                if sort([i1, i2, i3]) == unique_monomial(j,:)
                    polycoeff(j) = polycoeff(j) + T1(i1,i2,i3);
                end
            end
        end
    end
end

[sol, time] = SDDT(polycoeff, n, m);

