% Number of variables (indeterminates) in p
N = 5;
% m = 2;

x = msspoly('x',N);
prog = spotsosprog;
prog = prog.withIndeterminate(x);

g0 = 2*x(1) - x(2) + x(3) - 2*x(4) - 2*x(5);
g1 = (x(1) - 2)^2 - x(2)^2 - (x(3) - 1)^2 - (x(5) - 1)^2;
g2 = x(1)*x(3) - x(4)*x(5) + x(1)^2 - 1;
g3 = x(3) - x(2)^2 - x(4)^2 - 1;
g4 = x(1)*x(5) - x(2)*x(3) - 2;
g5 = 14 - x(1) - x(2) - x(3) - x(4) - x(5);
g6 = x(1);
g7 = x(2);
g8 = x(3);
g9 = x(4);
g10 = x(5);
g11 = x(1)^0;

g = {g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11};
% g = {g1, g11};
% g = {g1, g2, g3, g4, g5,g6, g7,  g8, g9,  0*x(1),  g11};

[prog,gamma] = prog.newFree(1);

f = g0 - gamma;

hier = 4;

[prog, f] = lasserre(prog, g, f, N, hier, x);
% [prog, f] = lasserre_dir(prog, g, f, N, hier, x);
% prog = prog.withSOS(f);
prog = prog.withPolyEqs(f);


% MOSEK options
% options = spot_sdp_default_options();
% options.solveroptions.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT'; % Use just the interior point algorithm to clean up
% options.solveroptions.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER'; % Don't use basis identification (it's slow)

%tic;
sol = prog.minimize(-gamma, @spot_mosek);
%time = toc;

% Optimal value
opt_sdsos = double(sol.eval(gamma));

disp(['Optimal value: ' num2str(opt_sdsos)])





function [prog, f] = lasserre(prog, g, f, N, hier, x)
    
    for i = 1: length(g)
        m = 2*floor(hier -deg(g{i})/2);
        vx =  monomials(x,0:m);
        [prog,sosvar] = prog.newFree(nchoosek(m+N, N));
        gtemp = sosvar'*vx;
        prog = prog.withSDSOS(gtemp);
%         if i == 10 || i == 9 || i== 8 || i ==7 || i ==6 || i ==7
%             f = f ;
%         else
%             f = f - gtemp*g{i};
%         end
        f = f - gtemp*g{i};
    end

end



