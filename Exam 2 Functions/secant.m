function xr=secant(f,x1,x2,tol,maxIt)
%xr=secant(f,x1,x2,tol,maxIt)
%
%Implements the Secant Method
%
%   f = function to find the root of
%   x1 = first x value
%   x2 = second x value
%   tol = stopping criteria tolerance
%   maxIt = maximum number of iterations
%   xr = found root
%
%ABH
%Spring 2023
if nargin()<3 || nargin()>5
error(['Inputs must include function and two starting points. Stopping ' ...
    'criteria and number of iterations are optional'])
elseif nargin()==3
tol=1e-10; maxIt=1e6;
else
maxIt=1e6;
end
y1=f(x1); y2=f(x2);
for ii=1:maxIt
xr=x2-y2*(x1-x2)/(y1-y2);
if abs(1-x2/xr)<tol
return
end
x1=x2; y1=y2;
x2=xr; y2=f(x2);
end
error('Reached maximum iterations and did not converge')
