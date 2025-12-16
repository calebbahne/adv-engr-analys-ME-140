function xr=modSecant(f,x0,delta,tol,maxIt)
%xr=modSecant(f,x0,delta,tol,maxIt)
%
%Implements the modified secant method
%
% f = function to find the root of
% x0 = initial value
% delta = value of delta
% tol = stopping criteria tolerance
% maxIt = maximum number of iterations
% xr = found root
%
%ABH
%Spring 2023
if nargin()<3 || nargin()>5
error(['Must provide function, initial x value and delta. Tolerance and ' ...
    'maximum iterations are optional.'])
elseif nargin()==3
tol=1e-10; maxIt=1e6;
else
maxIt=1e6;
end
for ii=1:maxIt
xr=x0-delta*x0*f(x0)/(f(x0+delta*x0)-f(x0));
if abs(1-x0/xr)<tol
return
end
x0=xr;
end
