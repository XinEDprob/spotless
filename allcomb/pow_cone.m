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
prob.cones.conepar= [0.5 0.5];
prob.cones.sub = [1 2 4 3 6 5];
prob.cones.subptr = [1 4];
[r,res]=mosekopt('maximize',prob);
% Display the primal solution.
res.sol.itr.xx'


