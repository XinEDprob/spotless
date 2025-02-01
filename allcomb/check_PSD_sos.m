% m = 2;
% n = 3;
% 
% x = msspoly('x',n); 
% p = x(1)^2 + 5*x(2)^2 + 3*x(3)^2;

m = 6;
n = 3;
x = msspoly('x',n); 
p = 10*x(1)^2*x(2)^2 + x(1)^6 + x(2)^4*x(3)^2 + x(2)^2*x(3)^4 - 1*x(1)^2*x(2)^2*x(3)^2;

% p =  1+ x(1)^6 + x(2)^6 + x(3)^6 - 3*x(1)^2*x(2)^2*x(3)^2;




%% DSOS program

% Initialize program
prog = spotsosprog;
prog = prog.withIndeterminate(x);

% New free variable gamma
[prog,gamma1] = prog.newFree(1);
[prog,gamma2] = prog.newFree(1);

% p = x(1)^6 + x(2)^6 + x(3)^6 - gamma1*x(1)^2*x(2)^2*x(3)^2 - (gamma2 )*x(1)^4*x(2)*x(3);

% p = x(1)^6 + x(2)^6 + x(3)^6 - gamma1*x(1)^2*x(2)^2*x(3)^2;

% p = x(1)^12 + x(2)^8*x(3)^4 + x(2)^4*x(3)^8 - 3*x(1)^4*x(2)^4*x(3)^4;
% p = x(1)^8 + x(2)^8 +  x(3)^8 - 1*x(1)^2*x(2)^3*x(3)^3;
% p = x(1)^2*x(2)^6 + x(2)^2*x(3)^6 +  x(1)^4*x(2)^4 - 1*x(1)^2*x(2)^2*x(3)^4;

% p = x(1)^6 + x(2)^4*x(3)^2 + x(2)^2*x(3)^4 - 1.4037*x(1)^2*x(2)^2*x(3)^2 - 1.4037*x(1)^4*x(2)*x(3);

% p = x(1)^6 + x(2)^6 + 0.99992*x(3)^6 - 1.4037*x(1)^2*x(2)^2*x(3)^2 - gamma2*x(1)^4*x(2)*x(3);


%%% SDP \neq SDSOS
% p = x(1)^12 + x(2)^12 + x(2)^12 +  10*x(1)^4*x(2)^8 - 1*x(1)^4*x(2)^4*x(3)^4 - gamma2*x(2)^6*x(3)^6 - gamma2*x(1)^8*x(2)^2*x(3)^2;;

% p = x(1)^12 + x(2)^12 + x(2)^12 +   10*x(2)^6*x(3)^6 - 1*x(1)^4*x(2)^4*x(3)^4 - gamma2*x(1)^8*x(2)^2*x(3)^2;
% 
% p = x(1)^12 + x(2)^8*x(3)^4 + x(2)^4*x(3)^8 +   10*x(2)^6*x(3)^6 - 1*x(1)^4*x(2)^4*x(3)^4 - gamma2*x(1)^8*x(2)^2*x(3)^2;


% p = x(1)^12 + x(2)^12 + x(2)^12  - 1*x(1)^4*x(2)^4*x(3)^4;

% p = x(1)^12 + x(2)^4*x(3)^8 + x(2)^8*x(3)^4  - x(1)^4*x(2)^4*x(3)^4;

% promising
% p =  1+  x(1)^6 + x(2)^6 + x(3)^6 + x(2)^2*x(3)^2 -gamma2*x(2)^3*x(1)^2    - 2.9*x(1)^2*x(2)^2*x(3)^2  ;


% p =   1+ x(1)^6 + x(2)^4*x(3)^2 + 10*x(2)^2*x(3)^4 -x(2)^2*x(3)^2 - gamma2*x(1)^2 - 1*x(1)^2*x(2)^2*x(3)^2 ;

% p =   1+ x(1)^6 + x(2)^4*x(3)^2 + x(2)^2*x(3)^4  - 3*x(1)^2*x(2)^2*x(3)^2 ;

% p =  1+ x(1)^6 + x(2)^4*x(3)^2 + x(2)^2*x(3)^4 + x(2)^2*x(3)^2 -0.62*x(1)^2*x(2)^2*x(3)    - 2.9*x(1)^2*x(2)^2*x(3)^2  ;

% p =   x(1)^6 + x(2)^4*x(3)^2 + x(2)^2*x(3)^4 + x(2)^2*x(3)^2 -0.034926*x(1)^2*x(2)^2*x(3)    - 2.9*x(1)^2*x(2)^2*x(3)^2  ;
% new
p =   x(1)^6 + x(2)^4*x(3)^2 + x(2)^2*x(3)^4 + x(2)^2*x(3)^2 + x(1)^2 -0.5*x(1)^2*x(2)^2*x(3)    - 2.9*x(1)^2*x(2)^2*x(3)^2  ;

%%%% check
% p =    x(1)^6 + x(2)^6 + x(3)^6 + x(2)^2*x(3)^2 -0.034926*x(2)^3*x(1)^2    - 2.9*x(1)^2*x(2)^2*x(3)^2  ;
% new
% p =    x(1)^6 + x(2)^6 + x(3)^6 + x(2)^2*x(3)^2 + x(1)^2 - 0.5*x(2)^3*x(1)^2    - 2.9*x(1)^2*x(2)^2*x(3)^2  ;



vx = monomials(x,0:m/2);

[prog,cp] = prog.newPSD(length(vx));

full_p = vx'*cp*vx;

prog = prog.withPolyEqs(p - full_p);

% prog = prog.withEqs(gamma1 - gamma2);


% [prog,ind1] = prog.withSDSOS(p);
% 


% MOSEK options
options = spot_sdp_default_options();
options.solveroptions.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT'; % Use just the interior point algorithm to clean up
options.solveroptions.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER'; % Don't use basis identification (it's slow)

% Solve program
%tic;
% sol = prog.minimize(gamma2 , @spot_mosek, options);
sol = prog.minimize(cp(3), @spot_mosek, options);
%time = toc;
% Optimal value
opt_dsos = double(sol.eval(gamma2));

disp(['Optimal value (DSOS): ' num2str(opt_dsos)])
Q = sol.eval(cp);
%disp(['Solution time: ' num2str(time)])

% if (sol.isPrimalFeasible) && (sol.isDualFeasible)
%     isdsos = true;
%     Q = sol.eval(cp);
%     phi = sol.gramMonomials{ind1};
% else
%     isdsos = false;
%     Q = [];
%     phi = [];
% end


%%
% syms x1 x2 x3
% syms z1 z2 z3
% z = [z1^3, z2^3, z3^3, z1*z2^2, z1^2*z2, z1*z3^2, z1^2*z3, z2*z3^2, z2^2*z3, z1^2*z2^2*z3^2];
% 
% x = subs(z, [z1,z3,z3], [x1,x2^(2/3)*x3^(1/3), x2^(1/3)*x3^(2/3)]);
% 
% A = round(10*rand(length(x), length(x)));
% 
% g = x*(A'*A)*x';