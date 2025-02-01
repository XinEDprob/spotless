% to test the relationship between DD induced polynomials and SDD induced
% polynomials


% Number of variables (indeterminates) in p
N = 5;
x = msspoly('x',N); % x of dimension N

%randn('state',4);


%% DSOS program

% Initialize program
prog = spotsosprog;
prog = prog.withIndeterminate(x);

% New free variable gamma
[prog,gamma] = prog.newFree(1);

% DSOS constraint
% prog = prog.withDSOS((p - gamma*(x'*x)^2)); % Only line that changes between DSOS,SDSOS,SOS programs
% prog = prog.withDSOS((p - gamma*(x(1)*x(1)*x(1)*x(1) + x(2)*x(2)*x(2)*x(2) + x(3)*x(3)*x(3)*x(3)+...
%     x(4)*x(4)*x(4)*x(4) + x(5)*x(5)*x(5)*x(5) ))); % Only line that changes between DSOS,SDSOS,SOS programs

%prog = prog.withDSOS(( 100 * (x(1)*x(1)*x(1)*x(2)  ) - gamma*(x(1)*x(1)*x(1)*x(1) + x(2)*x(2)*x(2)*x(2) + x(3)*x(3)*x(3)*x(3)+...
%     x(4)*x(4)*x(4)*x(4) + x(5)*x(5)*x(5)*x(5) )));

prog = prog.withDSOS( 100 *x(1)*x(3)*x(3)*x(3)*x(3)*x(3) - gamma*(x(1)*x(1)*x(1)*x(1)*x(1)*x(1) + x(2)*x(2)*x(2)*x(2)*x(2)*x(2)+ x(3)*x(3)*x(3)*x(3)*x(3)*x(3) ));


% MOSEK options
options = spot_sdp_default_options();
options.solveroptions.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT'; % Use just the interior point algorithm to clean up
options.solveroptions.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER'; % Don't use basis identification (it's slow)

% Solve program
sol = prog.minimize(-gamma, @spot_mosek, options);

% Optimal value
opt_dsos = double(sol.eval(gamma));

disp(['Optimal value (DSOS): ' num2str(opt_dsos)])