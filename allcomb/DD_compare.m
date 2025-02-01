%% to compare the optimization on DD matrices and DD tensors induced polynomials

%% Optimization on DD matrices induced polynomials

% Number of variables (indeterminates) in p
N = 20;
m  = 4;

x = msspoly('x',N); % x of dimension N

randn('state',10);

vx = monomials(x,m:m); % Degree 4 homogeneous
% Generate random polynomial
cp = randn(1,length(vx));
%cp = ones(1, length(vx));
p = cp*vx;

%% DSOS program

% Initialize program
prog = spotsosprog;
prog = prog.withIndeterminate(x);

% New free variable gamma
[prog,gamma] = prog.newFree(1);

% DSOS constraint
diagonal = 0;
for i = 1:N
    diagonal = diagonal + x(i)^m;
    
end



prog = prog.withDSOS(p - gamma * diagonal);

% prog = prog.withDSOS(p  - gamma);


% prog = prog.withDSOS(( 100 * (x(1)*x(1)*x(2)*x(2)*x(3)*x(3)  ) - gamma*(x(1)*x(1)*x(1)*x(1)*x(1)*x(1) + x(2)*x(2)*x(2)*x(2)*x(2)*x(2)+ x(3)*x(3)*x(3)*x(3)*x(3)*x(3) )));


% MOSEK options
options = spot_sdp_default_options();
options.solveroptions.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT'; % Use just the interior point algorithm to clean up
options.solveroptions.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER'; % Don't use basis identification (it's slow)

% Solve program
%tic;
sol = prog.minimize(-gamma, @spot_mosek, options);
%time = toc;
% Optimal value
opt_dsos = double(sol.eval(gamma));

disp(['Optimal value (DSOS): ' num2str(opt_dsos)])
%disp(['Solution time: ' num2str(time)])