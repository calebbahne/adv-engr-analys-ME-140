function xr=modSecant(f,x,delta,tol,maxIt)
%xr=modSecant(f,x,delta,tol,maxIt)
%
%Implements the modified secant method
%
%f is the function to find the root of
%x is the initial value
%delta is the value of delta
%tol is the stopping criteria tolerance
%maxIt is the maximum number of iterations
%xr is the found root
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
xr=x-delta*x*f(x)/(f(x+delta*x)-f(x));
if abs(1-x/xr)<tol
return
end
x=xr;
end
