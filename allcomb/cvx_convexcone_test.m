cvx_begin 
variables x1 x2 x3 x4 x5;
maximize(0)

subject to
x1*x2 >= power(x5, 2);
cvx_end
cvx_optval