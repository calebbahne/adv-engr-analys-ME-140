%% Caleb Bahne - ME 140 Exam 2
% 10.28.25
clc; clear;

%% Problem 21
% Bisection Method
f = @(x) -8-12*x+8*x.^2-3.25*x.^3;
fplot(f, [-1 0]) % looks to be between -1 and 0
xl = -1; xu = 0; es = 0.01;
[root,relErr,iter] = bisect(f,xl,xu,es)

% False Position Method
%   help falsepos
[root_fp,froot_fp,ea,iter_fp]=falsepos(f,xl,xu,es)

% Bisect Approx Error
f05 = f(-0.5)
fn1 = f(-1)
f1 = f(0)

% False pos approx err
ea   = 100;
xr   = xl;
xrold = xr;
fl = f(xl);
fu = f(xu);
xr = xu - fu*(xl - xu)/(fl - fu)   % false-position formula
ea = abs((xr - xrold)/xr) * 100
f_xr = f(xr)
xu = xr;
fu = f(xu);

xrold = xr;
xr = xu - fu*(xl - xu)/(fl - fu)
ea = abs((xr - xrold)/xr) * 100

%% Problem 22
% newtraph
func = @(x) x.^3+4*x.^2-10*x-5.2;
dfunc = @(x) 3*x.^2+8*x-10; 
fplot(func, [-7 3])
es = 0.1; xr = 3.5;
% help newtraph
[root_nr,ea,iter_nr]=newtraph(func,dfunc,xr,es)

% mod secant
help modSecant
delta = 0.001;
root_ms = modSecant(func,xr,delta,es)

% find all roots
coeff = [1 4 -10 -5.2];
roots(coeff)

%% Problem 23 - Golden Section Search
%    help goldmin
f = @(x) 2*x.^4+x.^3+10*x.^2+4*x;
g = @(x) -2*x.^4-x.^3-10*x.^2-4*x;
[x_min,fx,ea,iter]=goldmin(f,xl,xu,es)
[x_max,fx,ea,iter]=goldmin(g,xl,xu,es)

% help fminbnd
X = fminbnd(f, xl, xu) % min
X_max = fminbnd(g, xl, xu)

%% Problem 24 - Linear System
A = [2 -8 3; -4 1 7; 5 3 -2];
b = [40; -30; 50];
C = A\b

%% Problem 25 - Gauss Seidel
A = [2 -8 3; -4 1 7; 5 3 -2];
b = [20; -15; 25];
C = A\b

help GaussSeidel
[x_gs,ea,iter_gs] = GaussSeidel(A,b) % diverges
 
% diagonally dominant
A_dd = [5 3 -2; 2 -8 3; -4 1 7];
b_dd = [25; 20; -15];
[x_gsdd,ea,iter_gs] = GaussSeidel(A_dd,b_dd)

% with relaxation
% help GaussSeidelRelax
lambda = 0.01;
[x,ea,iter] = GaussSeidelRelax(A,b,lambda)