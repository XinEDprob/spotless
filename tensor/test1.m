prog = spotprog;
[prog,gamma] = prog.newPos(1);
prog = prog.withPos(gamma - 2);

% options = spot_sdp_default_options();
% options.solveroptions.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT'; % Use just the interior point algorithm to clean up
% options.solveroptions.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER'; % Don't use basis identification (it's slow)

sol = prog.minimize(gamma, @spot_mosek);
opt_sdsos = double(sol.eval(gamma));

disp(['Optimal value (DSOS): ' num2str(opt_sdsos)])
