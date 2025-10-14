function [x,fx,ea,iter]=goldmin(f,xl,xu,es,maxit,varargin)
% goldmin: minimization golden section search
% [x,fx,ea,iter]=goldmin(f,xl,xu,es,maxit,p1,p2,...):
% uses golden section search to find the minimum of f
% input:
%   f = name of function
%   xl, xu = lower and upper guesses
%   es = desired relative error (default = 0.0001%)
%   maxit = maximum allowable iterations (default = 50)
%   p1,p2,... = additional parameters used by f
% output:
%   x = location of minimum
%   fx = minimum function value
%   ea = approximate relative error (%)
%   iter = number of iterations
% Can find max w/ g(x) = -f(x) -> max of f(x) = -g(x)

if nargin<3,error('at least 3 input arguments required'),end
if nargin<4||isempty(es), es=0.0001;end
if nargin<5||isempty(maxit), maxit=50;end
phi=(1+sqrt(5))/2; d = (phi-1)*(xu - xl);
iter = 0; x1 = xl + d; x2 = xu - d;
f1 = f(x1,varargin{:}); f2 = f(x2,varargin{:});
while(1)
xint= xu - xl;
if f1 < f2
xopt = x1; xl = x2; x2 = x1; f2 = f1;
x1 = xl + (phi-1)*(xu-xl); f1 = f(x1,varargin{:});
else
xopt = x2; xu = x1; x1 = x2; f1 = f2;
x2 = xu - (phi-1)*(xu-xl); f2 = f(x2,varargin{:});
end
iter=iter+1;
if xopt~=0, ea = (2 - phi) * abs(xint / xopt) * 100;end
if ea <= es || iter >= maxit,break,end
end
x=xopt;fx=f(xopt,varargin{:});
