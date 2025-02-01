% Modified based on the code from sdsos_example.m 



% Number of variables (indeterminates) in p
N = 5  ;
m = 4;
x = msspoly('x',N); % x of dimension N

randn('state',10);

n_var = nchoosek(N+m-1, m);

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

Qt_inv = containers.Map;
for i = 1:DimT
    for i1 = 1:N
        for i2 = i1:N
            for i3 = i2:N
                for i4 = i3:N
                    if Qt(i1,i2,i3,i4) == i
                        Qt_inv(int2str(i)) = [i1,i2,i3,i4];
                    end
                end
            end
        end
    end
end

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
        prog = prog.withPos(cp(n_var + 1 - Qt(i,i,i,i)) - gamma -  ConsDD);
    end
    
    diagonal_term = ones(N,1);
    for i = 1:N
        diagonal_term(i) = Qt(i,i,i,i);
    end
    
    
    for i = 1:n_var
        if ~ismember(i, diagonal_term)
            index = Qt_inv(int2str(i));
            a = unique(index);
            out = [histc(index(:),a)];
            prog = prog.withPos(Auxt(i) - cp(n_var + 1 - i)/multinomial(4,out));
            prog = prog.withPos(Auxt(i) + cp(n_var +1 - i)/multinomial(4,out));
        end
    end
end


%% Solve

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
