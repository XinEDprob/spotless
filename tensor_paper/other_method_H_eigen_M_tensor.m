%%example 3.1
% A(1,:,:) = [7, 0, 0;
%     0, -0.5, -2;
%     0, -1, -2];
% A(2,:,:) = [-1, -5.8, -2;
%     0, 12, 0;
%     0, 0, -0.5];
% A(3,:,:) = [-1, -2, 0;
%     0, -1, -3;
%     0, -3, 50];

%%example 3.2
% A(1,:,:) = [12, 0, -0.3;
%     0, 0, -2;
%     0, -1, -1.5];
% A(2,:,:) = [-0.5, -4.8, -8;
%     0, 30, 0;
%     -1, 0, -0.5];
% A(3,:,:) = [0, -3, -1;
%     0, -1, -3.5;
%     -1, -3, 15];

%%symmetrilization of example 3.1
A(1,:,:) = [7, -1/3, -1/3;
    -1/3, -2.1, -7/6;
    -1/3, -7/6, -2/3];
A(2,:,:) = [-1/3, -2.1, -7/6;
    -2.1, 12, -1/3;
    -7/6, -1/3, -6.5/3];
A(3,:,:) = [-1/3, -7/6, -2/3;
    -7/6, -1/3, -6.5/3;
    -2/3, -6.5/3, 50];




% symmetrilization of example 3.2
% A(1,:,:) = [12, -0.5/3, -0.1;
%      -0.5/3, -1.6, -2.5;
%     -0.1, -2.5, -3.5/3];
% A(2,:,:) = [ -0.5/3, -1.6, -2.5;
%     -1.6, 30, -1/3;
%     -2.5, -1/3, -7/3];
% A(3,:,:) = [-0.1, -2.5, -3.5/3;
%     -2.5, -1, -7/3;
%     -3.5/3, -7/3, 15];

n = 3;
%% bounds 1
% lower bound
max_1 = -Inf;
min_1 = Inf;
for i = 1:n
    max_1 = max(max_1, sum(sum(A(i,:,:))));
    min_1 = min(min_1, sum(sum(A(i,:,:))));
end
fprintf('bounds 1\n')
fprintf('min: %d\n',min_1);
fprintf('max: %d\n',max_1);

%% bounds 2
ri = zeros(n,1);
rij = zeros(n);
for i = 1:n
    ri(i) = sum(sum(abs(A(i,:,:)))) - abs(A(i,i,i));
end


for i = 1:n
    for j = 1:n
        if i ~= j
            rij(i,j) = ri(i) - abs(A(i,j,j));
        else
            rij(i,j) = ri(i);
        end 
    end
end

Rij = zeros(n);
for i = 1:n
    for j = 1:n
        Rij(i,j) = (A(i,i,i) - A(j,j,j) - rij(i,j))^2 - 4*A(i,j,j)*ri(j);
    end
end


max_2 = -Inf;
min_2 = Inf;
for i = 1:n
    for j = 1:n
        if i~= j
            temp = A(i,i,i) + A(j,j,j) - rij(i,j) - sqrt(Rij(i,j));
            min_2 = min(min_2, temp);
            max_2 = max(max_2, temp);
        end
    end
end

fprintf('bounds 2\n')
fprintf('min: %d\n',min_2/2);
fprintf('max: %d\n',max_2/2);

%% bounds 3
Mij = zeros(n);
for i = 1:n
    for j = 1:n
        if i == j
            Mij(i,j) = ri(i);
        else
            Mij(i,j) = abs(A(i,j,j));
        end
    end
end

rM_i = zeros(n,1);
for i = 1:n
    rM_i(i) = sum(Mij(i,:)) - Mij(i,i);
end

r_tilde_i = zeros(n,1);
for i = 1:n
    r_tilde_i(i) = ri(i) - rM_i(i);
end

R_tilde_ij = zeros(n);
for i = 1:n
    for j = 1:n
        R_tilde_ij(i,j) = (A(i,i,i) - A(j,j,j) - r_tilde_i(i))^2 + 4*rM_i(i)*ri(j);
    end
end




max_3 = -Inf;
min_3 = Inf;
for i = 1:n
    for j = 1:n
        if i~= j
            temp = A(i,i,i) + A(j,j,j) - r_tilde_i(i) - sqrt(R_tilde_ij(i,j));
            min_3 = min(min_3, temp);
            max_3 = max(max_3, temp);
        end
    end
end

fprintf('bounds 3\n')
fprintf('min: %d\n',min_3/2);
fprintf('max: %d\n',max_3/2);

%% BOUNDS 4 

max_4 = -Inf;
min_4 = Inf;
for i = 1:n
    temp_max = -Inf;
    temp_min = Inf;
    for j = 1:n
        if i~= j
            temp = A(i,i,i) + A(j,j,j) - rij(i,j) - sqrt(Rij(i,j));
            temp_min = min(temp_min, temp);
            temp_max = max(temp_max, temp);
        end
    end
    max_4 = max(max_4, temp_min);
    min_4 = min(min_4, temp_max);
end

fprintf('bounds 4\n')
fprintf('min: %d\n',min_4/2);
fprintf('max: %d\n',max_4/2);


%% BOUNDS 5

r_theta_i = zeros(n,1);
for i = 1:n
    r_theta_i(i) = 0;
    for i2 = 1:n
        for i3 = 1:n
            if i2 == i || i3 == i
                r_theta_i(i) = r_theta_i(i) + abs(A(i,i2,i3));
            end
        end
    end
    r_theta_i(i) = r_theta_i(i) - abs(A(i,i,i));
end

r_theta_invi = zeros(n,1);
for i = 1:n
    r_theta_invi(i) = ri(i) - r_theta_i(i);
end

Delta_ij = zeros(n);
for i = 1:n
    for j= 1:n
        Delta_ij(i,j) = 0.5*(A(i,i,i) + A(j,j,j) - r_theta_i(i) - ...
            sqrt((A(i,i,i) - A(j,j,j) - r_theta_i(i))^2 + 4*r_theta_invi(i)*ri(j)));
    end
end


max_5 = -Inf;
min_5 = Inf;
for i = 1:n
    for j = 1:n
        if i~= j
            temp = Delta_ij(i,j);
            min_5 = min(min_5, temp);
            max_5 = max(max_5, temp);
        end
    end
end

fprintf('bounds 5\n')
fprintf('min: %d\n',min_5);
fprintf('max: %d\n',max_5);