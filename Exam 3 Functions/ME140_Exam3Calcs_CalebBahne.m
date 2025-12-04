%% ME 140 Exam 3 - Caleb Bahne
%% Problem 21

A = [-3 5; 110 -302];
% help eig
[V,D] = eig(A)

% polynomial method
p = poly(A) % char polynomial coeffs
% help poly
d = roots(p) % eigenvals
%% Problem 22

W = [65 72 76 80 83 85];
A = [2.06 2.12 2.14 2.18 2.20 2.25];

[a r2] = linregr(W,A)
% help linregr
% log10(A) = log10(alph) + beta*log10(W)
[a_pl, r2_pl] = linregr(log10(W), log10(A))
alpha = 10^(a_pl(2))
beta  = a_pl(1)

% check results
hold on;
fplot(@(x) alpha*x^beta, [65 85]); hold off;

% nonlinear method doublecheck
nlm = fitnlm(W,A, @(b,W) b(1)*W.^b(2), [1 1])
% seems consistent with transform and the Excel plot below
%% 
% 
%% Problem 23

x = [1 2 3 4 5 6 7];
y = [3.5 1.6 1.3 0.85 0.76 1.4 0.54];

% Newton interp poly
% help Newtint
xx = 1.7;
yint = Newtint(x,y,xx)

% Lagrange (no exclusion of the last point)
help Lagrange
xx = 1.5;
yint_L = Lagrange(x,y,xx)

% plot to evaluate the last value
plot(x,y,'*')

% Lagrange (exclude last point)
x1 = [1 2 3 4 5];
y1 = [3.5 1.6 1.3 0.85 0.76];
yint_L1 = Lagrange(x1,y1,xx)

% Relative error
re = abs(yint_L1-yint_L)/yint_L1

% polyint
polyint(x,y,xx) % similar to the Lagrange way
%% Problem 24

% Romberg integration
f = @(x) x.^2 + 4*x.^(-2) - 4;
a = 1; b = 4; es = 1;
% help romberg
[q,ea,iter] = romberg(f,a,b,es)

% integral
q1 = integral(f,a,b)