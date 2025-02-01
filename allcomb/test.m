clear prob;
[r, res] = mosekopt('symbcon');
% Specify the non-conic part of the problem.
prob.c = [-1 0 0 1 1 0];
prob.a = [1 1 0.5 0 0 0];
prob.blc = [2.0];
prob.buc = [2.0];
prob.blx = [-inf -inf -inf -inf -inf 1.0];
prob.bux = [ inf inf inf inf inf 1.0];
% Specify the cones.
prob.cones.type = [res.symbcon.MSK_CT_PPOW res.symbcon.MSK_CT_PPOW];
prob.cones.conepar= [0.2 0.4];
prob.cones.sub = [1 2 4 3 6 5];
prob.cones.subptr = [1 4];
[r,res]=mosekopt('maximize',prob);
% Display the primal solution.
res.sol.itr.xx'

%% 
clear;
prog = spotsosprog;

[prog,x] = prog.newPow(3,2);
%[prog, x] = prog.newRLor(3,2);
%[prog, x] = prog.newFree(3,2);

prog= prog.withPos(x(1,1) + x(2,1) + 0.5*x(1,2) - 2);
prog= prog.withPos(2 - x(1,1) - x(2,1) - 0.5*x(1,2) );

prog= prog.withPos(x(2,2) - 1);
prog= prog.withPos(1 -x(2,2));

prog.powpar = [0.2, 0.4];
options = spot_sdp_default_options();
options.solveroptions.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT'; % Use just the interior point algorithm to clean up
options.solveroptions.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER'; % Don't use basis identification (it's slow)

sol = prog.minimize( - x(3,1) -x(3,2) + x(1,1) , @spot_mosek, options);
%sol = prog.minimize( - sqrt(2)/2*x(3,1) -sqrt(2)/2*x(3,2) + x(1,1) , @spot_mosek, options);


% Optimal value
opt_sdsos = double(sol.eval(- x(3,1) -x(3,2) + x(1,1)));

disp(['Optimal value (DSOS): ' num2str(opt_sdsos)])