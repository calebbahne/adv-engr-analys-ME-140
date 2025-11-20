%% In-Class 10
% find the derivative of exp(-0.5*x)

clc; clear;
format long;

f = @(x) exp(-0.5*x);
help rombdiff
[d,ea,iter]=rombdiff(f,1);
q = d % derivative at x0 = 1

% evaluate using the diff builtin function
syms x
y = exp(-0.5*x);
df = diff(f, x);

% evaluate at x0 = 1
x0 = 1;
df_at_x0 = subs(df,x,x0);
df_eval = double(df_at_x0)