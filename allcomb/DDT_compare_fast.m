% Modified based on the code from sdsos_example.m 



% Number of variables (indeterminates) in p
N = 20  ;
m = 4;
x = msspoly('x',N); % x of dimension N

randn('state',10);

vx = monomials(x,m:m); % Degree 4 homogeneous
% Generate random polynomial
cp = randn(1,length(vx));
%cp = ones(1, length(vx));
p = cp*vx ;

%% Tensor experiment

% In this experiment, we will have the basis just from our problem. In the
% future, we may thik about writing a function to read and build the basis
% automatically from the problem

% Initialize program
prog = spotsosprog;
prog = prog.withIndeterminate(x);

% New free variable gamma
[prog,gamma] = prog.newFree(1);

% variables in the tensor matrix
DimT = length(vx);
[prog, Qtensor_vec] = prog.newFree(DimT);
%Qt = mss_v2sT(Qtensor_vec, N, m);
Qt = symtensor([1: DimT], m, N);
[prog, Auxt] = prog.newFree(DimT);
%expr = p - gamma*(x'*x)^2 ;
%expr = p - gamma*(x(1)*x(1)*x(1)*x(1) + x(2)*x(2)*x(2)*x(2) + x(3)*x(3)*x(3)*x(3)+ x(4)*x(4)*x(4)*x(4) + x(5)*x(5)*x(5)*x(5));
%expr = 100 * x(1)*x(1)*x(1)*x(3) + 100 * x(2)*x(2)*x(2)*x(1) - gamma*(x(1)*x(1)*x(1)*x(1) + x(2)*x(2)*x(2)*x(2) + x(3)*x(3)*x(3)*x(3)+ x(4)*x(4)*x(4)*x(4) + x(5)*x(5)*x(5)*x(5));
%exp = ( 100 * (x(1)*x(1)*x(2)*x(2)*x(3)*x(3)  ) - gamma*(x(1)*x(1)*x(1)*x(1)*x(1)*x(1) + x(2)*x(2)*x(2)*x(2)*x(2)*x(2)+ x(3)*x(3)*x(3)*x(3)*x(3)*x(3) ))

diagonal = 0;
for i = 1:N
    diagonal = diagonal + x(i)^m;
    
end

expr = p - gamma * diagonal; 

% expr = p + 15*diagonal - gamma; 

%% DD tensor constraints
% 4th order situation
if m == 4
    for i = 1:N
        ConsDD = 0;
        for i2 = 1:N
            for i3 = 1:N
                for i4 = 1:N
                    ConsDD = ConsDD + Auxt(Qt(i,i2,i3,i4));
                end
            end
        end
        ConsDD = ConsDD - Auxt(Qt(i,i,i,i));
        prog = prog.withPos(Qtensor_vec(Qt(i,i,i,i)) -  ConsDD);
    end


    for i = 1:length(Qtensor_vec)
        if ~ismember(i, [Qt(1,1,1,1),Qt(2,2,2,2),Qt(3,3,3,3),Qt(4,4,4,4),Qt(5,5,5,5)])
            prog = prog.withPos(Auxt(i) - Qtensor_vec(i));
            prog = prog.withPos(Auxt(i) + Qtensor_vec(i));
        end
    end

    % coefficient constraint
    Concoeff = 0;
    for i1 = 1:N
        for i2 = 1:N
            for i3 = 1:N
                for i4 = 1:N
                    Concoeff = Concoeff + Qtensor_vec(Qt(i1,i2,i3,i4)) * x(i1) * x(i2) * x(i3) * x(i4);
                end
            end
        end
    end
end

%% 6th order situation

if m == 6
    for i = 1:N
        ConsDD = 0;
        for i2 = 1:N
            for i3 = 1:N
                for i4 = 1:N
                    for i5 = 1:N
                        for i6 = 1:N
                            ConsDD = ConsDD + Auxt(Qt(i,i2,i3,i4,i5,i6));
                        end
                    end
                end
            end
        end
        ConsDD = ConsDD - Auxt(Qt(i,i,i,i,i,i));
        prog = prog.withPos(Qtensor_vec(Qt(i,i,i,i,i,i)) -  ConsDD);
    end


    diagonal_set = zeros(1,N);
    for i = 1:N
        diagonal_set(i) = Qt(i,i,i,i,i,i);
    end

    for i = 1:length(Qtensor_vec)
        if ~ismember(i, diagonal_set)
            prog = prog.withPos(Auxt(i) - Qtensor_vec(i));
            prog = prog.withPos(Auxt(i) + Qtensor_vec(i));
        end
    end

    % coefficient constraint
    Concoeff = 0;
    for i1 = 1:N
        for i2 = 1:N
            for i3 = 1:N
                for i4 = 1:N
                    for i5 = 1:N
                        for i6 = 1:N
                    Concoeff = Concoeff + Qtensor_vec(Qt(i1,i2,i3,i4,i5,i6)) *x(i1)*x(i2)*x(i3)*x(i4)*x(i5)*x(i6);
                        end
                    end
                end
            end
        end
    end
end


%% Solve

sosCnst = expr - Concoeff;

decvar = prog.variables;

A = diff(sosCnst,decvar);
b = subs(sosCnst,decvar,0*decvar);
[var,pow,Coeff] = decomp([b A].');

[prog,y] = prog.withEqs(Coeff'*[1;decvar]);
basis = recomp(var,pow,speye(size(pow,1))); 

%prog.withEqs(gamma - 2);

% MOSEK options
options = spot_sdp_default_options();
options.solveroptions.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT'; % Use just the interior point algorithm to clean up
options.solveroptions.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER'; % Don't use basis identification (it's slow)

%tic;
sol = prog.minimize(-gamma, @spot_mosek, options);
%time = toc;

% Optimal value
opt_sdsos = double(sol.eval(gamma));

disp(['Optimal value (DSOS): ' num2str(opt_sdsos)])
%disp(['Solution time: ' num2str(time)])
