function [d,ea,iter]=rombdiff(func,x,es,maxit,varargin)
% romberg: Romberg integration quadrature
% q = romberg(func,a,b,es,maxit,varargin):
% [d,ea,iter]=rombdiff(func,x,es,maxit,varargin)
% Romberg differentiation.
% input:
% func = name of function to be integrated
% es = desired relative error (default = 0.000001%)
% maxit = maximum allowable iterations (default = 50)
% p1,p2,... = additional parameters used by func
% output:
% q = integral estimate
% ea = approximate relative error (%)
% iter = number of iterations
if nargin<2,error('at least 2 input arguments required'),end
if nargin<3|isempty(es), es=0.000001;end
if nargin<4|isempty(maxit), maxit=50;end
n = 1;
DY(1,1) = dydxnew(func,x,n,varargin{:});
iter = 0;
ea=100;
while iter<maxit
iter = iter+1;
n = 2^iter;
DY(iter+1,1) = dydxnew(func,x,n,varargin{:});
for k = 2:iter+1
j = 2+iter-k;
DY(j,k) = (4^(k-1)*DY(j+1,k-1)-DY(j,k-1))/(4^(k-1)-1);
end
if DY(1,iter+1)~=0,ea = abs((DY(1,iter+1)-DY(2,iter))/DY(1,iter+1))*100;end
if ea<=es, break; end
end
d = DY(1,iter+1);
