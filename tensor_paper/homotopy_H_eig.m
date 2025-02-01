%% Table 3.2 (methods from https://users.math.msu.edu/users/chenlipi/teneig.html)

% The followings are two examples, to test one of them, you just need to
% comment out another one
% symmetrized Example 3.1 from https://www.sciencedirect.com/science/article/pii/S002437951830404X
% A(1,:,:) = [7, -1/3, -1/3;
%     -1/3, -2.1, -7/6;
%     -1/3, -7/6, -2/3];
% A(2,:,:) = [-1/3, -2.1, -7/6;
%     -2.1, 12, -1/3;
%     -7/6, -1/3, -6.5/3];
% A(3,:,:) = [-1/3, -7/6, -2/3;
%     -7/6, -1/3, -6.5/3;
%     -2/3, -6.5/3, 50];


% symmetrized Example 3.2 from https://www.sciencedirect.com/science/article/pii/S002437951830404X
A(1,:,:) = [12, -0.5/3, -0.1;
     -0.5/3, -1.6, -2.5;
    -0.1, -2.5, -3.5/3];
A(2,:,:) = [ -0.5/3, -1.6, -2.5;
    -1.6, 30, -1/3;
    -2.5, -1/3, -7/3];
A(3,:,:) = [-0.1, -2.5, -3.5/3;
    -2.5, -1, -7/3;
    -3.5/3, -7/3, 15];


% ajacent tensor of 4-node 3-uniform hypergraph 
% A(1,:,:) = [0, 0, 0, 0;
%             0, 0, -1/2, 0;
%             0, -1/2, 0, -1/2;
%             0, 0, -1/2, 0];
% A(2,:,:) = [0, 0, -1/2, 0;
%             0, 0, 0, 0;
%             -1/2, 0, 0, 0;
%             0, 0, 0, 0];
% A(3,:,:) = [0, -1/2, 0, -1/2;
%             -1/2, 0, 0, 0;
%             0, 0, 0, 0;
%             -1/2, 0, 0, 0];
% A(4,:,:) = [0, 0, -1/2, 0;
%             0, 0, 0, 0;
%             -1/2, 0, 0, 0;
%             0, 0, 0, 0];
  

tic;
[lambda, V]  = heig(A,'symmetric');
toc
lambda