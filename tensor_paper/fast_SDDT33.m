 
function [opt_sdsos, time] = fast_SDDT33(cp)
 tic
 load SDDT_N3m3.mat gamma n_var diagonal_term Qt_inv prog Diag Offdiag

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

%% Solve

% MOSEK options
options = spot_sdp_default_options();
options.solveroptions.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT'; % Use just the interior point algorithm to clean up
options.solveroptions.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER'; % Don't use basis identification (it's slow)

% tic
sol = prog.minimize(-gamma, @spot_mosek, options);

time = toc;
disp(['Solution time: ' num2str(time)])
% Optimal value
opt_sdsos = double(sol.eval(gamma));

disp(['Optimal value (DSOS): ' num2str(opt_sdsos)])
%disp(['Solution time: ' num2str(time)])

end