x = [1 2 3 4 5 6 7];
f = [3.6 1.8 1.2 0.9 0.72 1.5 0.5149];
xx = 1.7;
% help Newtint
% help Lagrange

yint_N = Newtint(x,f,xx)
yint_L = Lagrange(x,f,xx)
% results are the same for both methods
