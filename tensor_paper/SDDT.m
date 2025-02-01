% Modified based on the code from sdsos_example.m 
% format long %load 
%% Test for minimum eigenvalue of M-tensor

% cp = [50, -6.5, -1, 12, -2, -7, -6.3, -1, -1, 7];

%  cp = [7, -1, -1, -6.3, -7, -2, 12, -1, -6.5, 50];

% cp = [15, -7, -1, 30, -3.5, -15,  -4.8, -0.3, -0.5, 12];

%  cp = [12, -0.5, -0.3, -4.8, -15, -3.5, 30, -1, -7, 15];

%  N = 3;
%  m = 3;
% hypergraph example %%%%%%%%%%%%%%
% cp = zeros(1, 20);
% cp(6) = -3;
% cp(9) = -3;
% N = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% test problem
% cp = [7, -1, -1, 0, 0, 0, 7, 0, 0, 7];

%example 1 for system of multilinear equations
% N = 3;
% m = 3;
% cp = polycoeff;

% example 2 for system of multilinear equations
% N = 3;
% m = 3;
% cp = polycoeff;

function[opt_sdsos, time] = SDDT(cp, N, m)

tic;
%% Number of variables (indeterminates) in p
% N = 4;
x = msspoly('x',N);
randn('state',10);


n_var = nchoosek(N+m-1, m);
%% Tensor experiment

% In this experiment, we will have the basis just from our problem. In the
% future, we may thik about writing a function to read and build the basis
% automatically from the problem

% Initialize program
prog = spotsosprog;
prog = prog.withIndeterminate(x);

% New free variable gamma
[prog,gamma] = prog.newFree(1);

%% SDD tensor constraints
% 4th order situation
if m == 4
    Qt = symtensor([1: n_var], m, N);

    unique_monomial = [];

    Qt_inv = containers.Map;

    i = 1;
    for i1 = 1:N
        for i2 = i1:N
            for i3 = i2:N
                for i4 = i3:N
                    Qt_inv(int2str(i)) = [i1,i2,i3,i4];
                    unique_monomial = [unique_monomial; [i1,i2,i3,i4]];
                    i = i+1;
                end
            end
        end
    end


    % construct the variable polynomial

    num_var_pow = 2 * 3 * nchoosek(N, 2)  + 3 * 3 * nchoosek(N , 3) + 4 *nchoosek(N, 4);
    [prog, diag_var] = prog.newPos(N, num_var_pow/N);

    num_pow = nchoosek(N, 2) + 2 * 2 * nchoosek(N, 2) + 2 * 3 * nchoosek(N , 3) + 3 *nchoosek(N, 4);
    %num_pow = nchoosek(N, 2) + 2 * 2 * nchoosek(N, 2) + 2 * 3 * nchoosek(N , 3)  ;
    [prog, pow_var] =  prog.newPow(3,num_pow);

    num_unique_monomial = nchoosek(N+m - 1, m);
    [prog, offdiag_var] = prog.newFree(num_unique_monomial);



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


    diagonal_term = ones(N,1);
    for j = 1:N
        diagonal_term(j) = Qt(j,j,j,j);
    end

    for j = 1:n_var
       if ~ismember(j, diagonal_term)
            index = Qt_inv(int2str(j));
            a = unique(index);
            out = [histc(index(:),a)];
            prog = prog.withEqs(offdiag_var(j) - cp(n_var + 1 - j));
    %         prog = prog.withPos(offdiag_var(i) + cp(n_var + 1 - i)/multinomial(4,out));
       end
        if ismember(j, diagonal_term)
            index = Qt_inv(int2str(j));
            a = unique(index);
            prog = prog.withPos(cp(n_var + 1 - j) - gamma - sum(diag_var(a,:)));
        end
    end

end
%% SDD tensor constraints
% 3th order situation
if m == 3
    Qt = symtensor([1: n_var], m, N);

    unique_monomial = [];

    Qt_inv = containers.Map;

    i = 1;
    for i1 = 1:N
        for i2 = i1:N
            for i3 = i2:N
                Qt_inv(int2str(i)) = [i1,i2,i3];
                unique_monomial = [unique_monomial; [i1,i2,i3]];
                i = i+1;
            end
        end
    end



    % construct the variable polynomial
    if N >= 3
        num_var_pow = 2 * 2 * nchoosek(N, 2)  + 3 *  nchoosek(N , 3);
    else
        num_var_pow = 2 * 2 * nchoosek(N, 2);
    end
    [prog, diag_var] = prog.newPos(N, num_var_pow/N);
    
    if N >=3
        num_pow =  2 *  nchoosek(N, 2) + 2 *  nchoosek(N , 3);
    else
        num_pow =  2 *  nchoosek(N, 2);
    end
    %num_pow = nchoosek(N, 2) + 2 * 2 * nchoosek(N, 2) + 2 * 3 * nchoosek(N , 3)  ;
    [prog, pow_var] =  prog.newPow(3,num_pow);

    num_unique_monomial = nchoosek(N+m - 1, m);
%     [prog, offdiag_var] = prog.newPos(num_unique_monomial);
    [prog, offdiag_var] = prog.newFree(num_unique_monomial);


    % hashtable for diag terms
    Diag = diag_var;
    
    % for off-diag terms
    Offdiag = offdiag_var;
    
    % record the num of diagonal variables used
    record_diag = ones(1,N);
    % record the num of pow constraints used 
    record_pow = 1;
    prog.powpar = [];
    for i  = 1: size(unique_monomial, 1)
        [C,ia,ic] = unique(unique_monomial(i,:));
        a_counts = accumarray(ic,1);
        if length(C) == 1
            prog = prog.withPos(offdiag_var(i));
        end

        if length(C) == 2
            if isequal(a_counts, [1;2]) || isequal(a_counts, [2;1])
                % index for the cubic term
                index3 = C(find(a_counts == 2 , 1 , 'first'));
                index1 = C(find(a_counts == 1 , 1 , 'first'));

                Diag(index3, record_diag(index3)) = pow_var(2, record_pow);
                Diag(index1, record_diag(index1)) = pow_var(1, record_pow);
                
                Offdiag(i) = pow_var(3,record_pow)/(nthroot(4,3)/3); 
%                 prog = prog.withPos(diag_var(index3, record_diag(index3)) - pow_var(2, record_pow));
%                 prog = prog.withPos(pow_var(2, record_pow) - diag_var(index3, record_diag(index3)));
%                 
%                 prog = prog.withPos(diag_var(index1, record_diag(index1)) - pow_var(1, record_pow));
%                 prog = prog.withPos(pow_var(1, record_pow) - diag_var(index1, record_diag(index1)));

%                 prog = prog.withPos(offdiag_var(i)*nthroot(4,3)/3 - pow_var(3,record_pow) );
%                 prog = prog.withPos(pow_var(3,record_pow) - offdiag_var(i)*nthroot(4,3)/3);


                prog.powpar = [prog.powpar 1/3];

                record_diag(index3) = record_diag(index3) + 1;
                record_diag(index1) = record_diag(index1) + 1;
                record_pow = record_pow + 1;
%                 prog = prog.withPos(-offdiag_var(i));
            end
        end
        
        % for case (1,2,3)
        if length(C) == 3
            index1 = C(1);
            index2 = C(2);
            index3 = C(3);
            
            Diag(index1, record_diag(index1)) = pow_var(1, record_pow);
            Diag(index2, record_diag(index2)) = pow_var(2, record_pow);
            
            
%             prog = prog.withPos(diag_var(index1, record_diag(index1)) - pow_var(1, record_pow));
%             prog = prog.withPos(pow_var(1, record_pow) - diag_var(index1, record_diag(index1)));
%             prog = prog.withPos(diag_var(index2, record_diag(index2)) - pow_var(2, record_pow));
%             prog = prog.withPos(pow_var(2, record_pow) - diag_var(index2, record_diag(index2)));

            prog = prog.withPos(pow_var(1, record_pow+1) - pow_var(3, record_pow));
            prog = prog.withPos(pow_var(3, record_pow) - pow_var(1, record_pow+1));       


            prog.powpar = [prog.powpar 0.5];

            
            Diag(index3, record_diag(index3)) = pow_var(2, record_pow+1);
            
            Offdiag(i) = 3*pow_var(3,record_pow + 1);
%             prog = prog.withPos(diag_var(index3, record_diag(index3)) - pow_var(2, record_pow + 1));
%             prog = prog.withPos(pow_var(2, record_pow + 1) - diag_var(index3, record_diag(index3)));

%             prog = prog.withPos(offdiag_var(i)/3  - pow_var(3,record_pow + 1) );
%             prog = prog.withPos(pow_var(3,record_pow + 1) - offdiag_var(i)/3 ) ;

            prog.powpar = [prog.powpar 2/3];

            record_diag(index1) = record_diag(index1) + 1;
            record_diag(index2) = record_diag(index2) + 1;
            record_diag(index3) = record_diag(index3) + 1;
            record_pow = record_pow + 2;
%             prog = prog.withPos(-offdiag_var(i));
        end


    end


    diagonal_term = ones(N,1);
    for j = 1:N
        diagonal_term(j) = Qt(j,j,j);
    end


    
    for j = 1:n_var
       if ~ismember(j, diagonal_term)
            index = Qt_inv(int2str(j));
            a = unique(index);
            out = [histc(index(:),a)];
            prog = prog.withEqs(Offdiag(j) - cp(j));
    %         prog = prog.withPos(offdiag_var(i) + cp(n_var + 1 - i)/multinomial(4,out));
       end
        if ismember(j, diagonal_term)
            index = Qt_inv(int2str(j));
            a = unique(index);
            prog = prog.withPos(cp(j) - gamma - sum(Diag(a,:)) - Offdiag(j));
        end
    end

end




%% Solve



% MOSEK options
options = spot_sdp_default_options();
options.solveroptions.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT'; % Use just the interior point algorithm to clean up
options.solveroptions.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER'; % Don't use basis identification (it's slow)

% disp('what?')
% tic
sol = prog.minimize(-gamma, @spot_mosek, options);

time = toc;
disp(['Solution time: ' num2str(time)])
% Optimal value
opt_sdsos = double(sol.eval(gamma));

disp(['Optimal value (DSOS): ' num2str(opt_sdsos)])
%disp(['Solution time: ' num2str(time)])

end
