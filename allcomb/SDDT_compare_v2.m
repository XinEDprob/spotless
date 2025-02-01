% Modified based on the code from sdsos_example.m 



% Number of variables (indeterminates) in p
N = 5;
m = 4;
x = msspoly('x',N); % x of dimension N

randn('state',10);

vx = monomials(x,m:m); % Degree 4 homogeneous
% Generate random polynomial
cp = randn(1,length(vx));
%cp = ones(1, length(vx));
p = cp*vx;

%p = 100*x(1) * x(1) *x(1) *x(1);
%% Tensor experiment

% In this experiment, we will have the basis just from our problem. In the
% future, we may thik about writing a function to read and build the basis
% automatically from the problem

% Initialize program
prog = spotsosprog;
prog = prog.withIndeterminate(x);

% New free variable gamma
[prog,gamma] = prog.newFree(1);

% variables in the tensor matrix;
%expr = p - gamma*(x'*x)^2;
%expr = p - gamma*(x(1)*x(1)*x(1)*x(1) + x(2)*x(2)*x(2)*x(2) + x(3)*x(3)*x(3)*x(3)+ x(4)*x(4)*x(4)*x(4) + x(5)*x(5)*x(5)*x(5));
%expr = 100 * x(1)*x(1)*x(1)*x(3) + 100 * x(2)*x(2)*x(2)*x(1) - gamma*(x(1)*x(1)*x(1)*x(1) + x(2)*x(2)*x(2)*x(2) + x(3)*x(3)*x(3)*x(3)+ x(4)*x(4)*x(4)*x(4) + x(5)*x(5)*x(5)*x(5));
%exp = ( 100 * (x(1)*x(1)*x(2)*x(2)*x(3)*x(3)  ) - gamma*(x(1)*x(1)*x(1)*x(1)*x(1)*x(1) + x(2)*x(2)*x(2)*x(2)*x(2)*x(2)+ x(3)*x(3)*x(3)*x(3)*x(3)*x(3) ))

diagonal = 0;
for i = 1:N
    diagonal = diagonal + x(i)^m;
    
end

expr = p - gamma * diagonal;



%% DD tensor constraints
% 4th order situation
if m == 4
    all_monomial = allcomb(1:N, 1:N, 1:N, 1:N);
end

unique_monomial = [1,1,1,1];

for i = 1: size(all_monomial,  1)
    indicator = 1;
    for j = 1:size(unique_monomial, 1)
        if ismember(all_monomial(i,:), perms(unique_monomial(j,:)), 'row')
            indicator = 0;
            break;
        end
    end
    if indicator == 1
        unique_monomial = [unique_monomial; all_monomial(i,:)];
    end
end


% construct the variable polynomial

num_var_pow = 2 * 3 * nchoosek(N, 2)  + 3 * 3 * nchoosek(N , 3) + 4 *nchoosek(N, 4);
[prog, diag_var] = prog.newPos(N, num_var_pow/N);

num_pow = nchoosek(N, 2) + 2 * 2 * nchoosek(N, 2) + 2 * 3 * nchoosek(N , 3) + 3 *nchoosek(N, 4);
%num_pow = nchoosek(N, 2) + 2 * 2 * nchoosek(N, 2) + 2 * 3 * nchoosek(N , 3)  ;
[prog, pow_var] =  prog.newPow(3,num_pow);

num_unique_monomial = (nchoosek(m + N, N) - nchoosek(m + N - 1, m - 1));
[prog, offdiag_var] = prog.newFree(num_unique_monomial);

% % hashtable for diag terms
% Diag = diag_var;
% 
% % for off-diag terms
% Offdiag = offdiag_var;

coeff = 0;
for i = 1:N
    coeff = coeff + sum(diag_var(i,:)) * x(i)^m;
end

for i = 1: size(unique_monomial, 1)
    coeff = coeff +  offdiag_var(i) * x(unique_monomial(i,1))* x(unique_monomial(i,2))* x(unique_monomial(i,3))* x(unique_monomial(i,4));

end


% record the num of diagonal variables used
record_diag = ones(1,N);
% record the num of socp variables used 
record_pow = 1;
prog.powpar = [];
for i  = 1: size(unique_monomial, 1)
    [C,ia,ic] = unique(unique_monomial(i,:));
    a_counts = accumarray(ic,1);
    if length(C) == 1
        prog = prog.withPos(offdiag_var(i));
    end
    
    if length(C) == 2
        if isequal(a_counts, [2;2])
            prog = prog.withPos(diag_var(C(1), record_diag(C(1))) - pow_var(1,record_pow));
            prog = prog.withPos(pow_var(1,record_pow) - diag_var(C(1), record_diag(C(1))));
            prog = prog.withPos(diag_var(C(2), record_diag(C(2))) - pow_var(2,record_pow));
            prog = prog.withPos(pow_var(2,record_pow) - diag_var(C(2), record_diag(C(2))));
            
            prog = prog.withPos(offdiag_var(i)/2 - pow_var(3,record_pow) );
            prog = prog.withPos(pow_var(3,record_pow) - offdiag_var(i)/2);
            
            record_diag(C(1)) = record_diag(C(1)) + 1;
            record_diag(C(2)) = record_diag(C(2)) + 1;
            record_pow = record_pow + 1;
            
            prog.powpar = [prog.powpar 0.5];
        end
        
        if isequal(a_counts, [1;3]) || isequal(a_counts, [3;1])
            % index for the cubic term
            index3 = C(find(a_counts == 3 , 1 , 'first'));
            index1 = C(find(a_counts == 1 , 1 , 'first'));
            
            prog = prog.withPos(pow_var(2, record_pow) - pow_var(1, record_pow));
            
            
            prog = prog.withPos(diag_var(index3, record_diag(index3)) - pow_var(2, record_pow));
            prog = prog.withPos(pow_var(2, record_pow) - diag_var(index3, record_diag(index3)));
            
            prog = prog.withPos(pow_var(1, record_pow+1) - pow_var(3, record_pow));
            prog = prog.withPos(pow_var(3, record_pow) - pow_var(1, record_pow+1));
            
            
            prog.powpar = [prog.powpar 2/3];
            
         
            prog = prog.withPos(diag_var(index1, record_diag(index1)) - pow_var(2, record_pow + 1));
            prog = prog.withPos(pow_var(2, record_pow + 1) - diag_var(index1, record_diag(index1)));
            
            prog = prog.withPos(offdiag_var(i)*nthroot(27,4)/4 - pow_var(3,record_pow + 1) );
            prog = prog.withPos(pow_var(3,record_pow + 1) - offdiag_var(i)*nthroot(27,4)/4 ) ;
            
            prog.powpar = [prog.powpar 3/4];

            
            record_diag(index3) = record_diag(index3) + 1;
            record_diag(index1) = record_diag(index1) + 1;
            record_pow = record_pow + 2;
            
        end
    end
    
    if length(C) == 3
        index2 = C(find(a_counts == 2 , 1 , 'first'));
        temp = find(a_counts == 1);
        index1 = C(temp(1));
        index11 = C(temp(2));
        
        prog = prog.withPos(diag_var(index2, record_diag(index2)) - pow_var(1, record_pow));
        
        prog = prog.withPos(diag_var(index1, record_diag(index1)) - pow_var(2, record_pow));
        prog = prog.withPos(pow_var(2, record_pow) - diag_var(index1, record_diag(index1)));
        
        prog = prog.withPos(pow_var(1, record_pow+1) - pow_var(3, record_pow));
        prog = prog.withPos(pow_var(3, record_pow) - pow_var(1, record_pow+1));       
        
        
        prog.powpar = [prog.powpar 2/3];
        
        prog = prog.withPos(diag_var(index11, record_diag(index11)) - pow_var(2, record_pow + 1));
        prog = prog.withPos(pow_var(2, record_pow + 1) - diag_var(index11, record_diag(index11)));      
        
        prog = prog.withPos(offdiag_var(i)*nthroot(3*3*6*6, 4)/12 - pow_var(3,record_pow + 1) );
        prog = prog.withPos(pow_var(3,record_pow + 1) - offdiag_var(i)*nthroot(3*3*6*6, 4)/12 ) ;
        
        prog.powpar = [prog.powpar 3/4];
        
        record_diag(index2) = record_diag(index2) + 1;
        record_diag(index11) = record_diag(index11) + 1;
        record_diag(index1) = record_diag(index1) + 1;
        record_pow = record_pow + 2;
        
    end
    if length(C) == 4
        index1 = C(1);
        index2 = C(2);
        index3 = C(3);
        index4 = C(4);

        prog = prog.withPos(diag_var(index1, record_diag(index1)) - pow_var(1, record_pow));
        prog = prog.withPos(pow_var(1, record_pow) - diag_var(index1, record_diag(index1)));
        prog = prog.withPos(diag_var(index2, record_diag(index2)) - pow_var(2, record_pow));
        prog = prog.withPos(pow_var(2, record_pow) - diag_var(index2, record_diag(index2)));
        
        prog = prog.withPos(pow_var(1, record_pow+1) - pow_var(3, record_pow));
        prog = prog.withPos(pow_var(3, record_pow) - pow_var(1, record_pow+1));       
        
        
        prog.powpar = [prog.powpar 0.5];
        
        
        prog = prog.withPos(diag_var(index3, record_diag(index3)) - pow_var(2, record_pow + 1));
        prog = prog.withPos(pow_var(2, record_pow + 1) - diag_var(index3, record_diag(index3)));
        
        prog = prog.withPos(pow_var(1, record_pow + 2) - pow_var(3, record_pow + 1));
        prog = prog.withPos(pow_var(3, record_pow + 1) - pow_var(1, record_pow + 2));
        
        prog.powpar = [prog.powpar 2/3];
        
        prog = prog.withPos(diag_var(index4, record_diag(index4)) - pow_var(2, record_pow + 2));
        prog = prog.withPos(pow_var(2, record_pow + 2) - diag_var(index4, record_diag(index4)));
        
        
        
        prog = prog.withPos(offdiag_var(i)/4  - pow_var(3,record_pow + 2) );
        prog = prog.withPos(pow_var(3,record_pow + 2) - offdiag_var(i)/4 ) ;

        prog.powpar = [prog.powpar 3/4];
        
        record_diag(index1) = record_diag(index1) + 1;
        record_diag(index2) = record_diag(index2) + 1;
        record_diag(index3) = record_diag(index3) + 1;
        record_diag(index4) = record_diag(index4) + 1;
        record_pow = record_pow + 3;
        
    end
    
end



%% 6th order situation






sosCnst = expr - coeff;

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

% Optimal value
opt_sdsos = double(sol.eval(gamma));

disp(['Optimal value (DSOS): ' num2str(opt_sdsos)])
%disp(['Solution time: ' num2str(time)])
