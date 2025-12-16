%% ME 140 Final Exam - Caleb Bahne
%% Problem 41

m = 70; t = 2; g = 9.81; cd = 0.3;
v = sqrt(g*m/cd)*tanh(sqrt(g*cd/m)*t)
%% Problem 42

% find roots between xl and xu
f = @(x) -6-10*x+9*x.^2-7*x.^3;
xl = -1; xu = 0; es = 0.01;
help bisect

[root,relErr,iter] = bisect(f,xl,xu,es)

fplot(f, [-2 2]); grid on;

% Find all roots
coeff = [-7 9 -10 -6];
all_roots = roots(coeff)
%% Problem 43

% Solve for x1, x2, x3
A = [3 -8 2; -4 1 9; 5 3 -4];
B = [25; -35; 66];
C = A\B

% Check ans
A(1,1)*C(1)+A(1,2)*C(2)+A(1,3)*C(3) % should be 25
A(2,1)*C(1)+A(2,2)*C(2)+A(2,3)*C(3) % -35
A(3,1)*C(1)+A(3,2)*C(2)+A(3,3)*C(3) % 66

% Use Gauss Seidel
help GaussSeidel
es = 0.5;
[x,ea,iter] = GaussSeidel(A,B,es) % diverges

% Get things diagonally dominant
A2 = [5 3 -4; 3 -8 2; -4 1 9];
B2 = [66; 25; -35];
[x,ea,iter] = GaussSeidel(A2,B2,es)
%% Problem 44

W = [50 62 71 80 83 85];
A = [1.85 1.94 2 2.1 2.13 2.15];

% Linearize for A = aW^b
% log10(A) = log10(alpha) + beta*log10(W)
[a, r2] = linregr(log10(W), log10(A))
alpha = 10^(a(2))
beta  = a(1)

% Check by plotting
f = alpha*W.^beta;
plot(W,A,'b',W,f,'g')

% Power law equation
nlm = fitnlm(W,A,@(b,W)b(1)*W.^b(2), [1 1])
%% Problem 45

% eigenvalues and vectors
help eig
A = [-4 6; 66 -313];
[V,D] = eig(A)

% Check with characteristic eq
P = poly(A)
lamda = roots(P)
%% Problem 46

% Derivative
syms x
f = 3*exp(-0.9*x);
df = diff(f)
dfunc = matlabFunction(df);
dfunc_2 = dfunc(2)

% Integral
f = @(x) (x-5./x).^4;
F = integral(f,1,3)

% check
hold off;
fplot(f,[1 3])
%% Extra Credit

tspan = [0 2];
y0 = 1; h = 0.001;
[t,y] = rk4(@dydtsys,tspan,y0,h);
plot(t,y)

% integrate
trapuneq(t,y)

% Notes
%   - you could also find a nonlinear model for y
%   - you could differentiate y at any point on [0 2]
%   - you could find the minimum of y on the interval