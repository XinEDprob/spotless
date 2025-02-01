N = 2;
x = msspoly('x',N);
prog = spotsosprog;
prog = prog.withIndeterminate(x);

g0 = x(1)^2* (4 - 2.1*x(1)^2 + 1/3*x(1)^4) + x(1)*x(2) +x(2)^2*(-4 + 4*x(2)^2);  
% g0 = x(1)^2 + x(2)^2  - 1;
g1 = 3 - x(1)^2 - x(2)^2;
g2 = x(1)^0;


g = {g1, g2};
% g = {g1, g11};


[prog,gamma] = prog.newFree(1);

f = g0 - gamma;

hier = 5;

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
        prog = prog.withSOS(gtemp);
        f = f - gtemp*g{i};
    end

end

function [prog, f] = lasserre_dir(prog, g, f, N, hier, x)
    
    for i = 1: length(g)
        m = 2*floor(hier -deg(g{i})/2);
        vx =  monomials(x,0:m);
        [prog,sospoly] = prog.newSOSPoly(vx);
%         gtemp = sosvar'*vx;
%         prog = prog.withSOS(gtemp);
        f = f - sospoly*g{i};
    end

end


