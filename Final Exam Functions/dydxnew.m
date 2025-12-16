function d = dydxnew(func,x,n,varargin)
a = x - x/n;
b = x + x/n;
d=(func(b,varargin{:})-func(a,varargin{:}))/(b-a);
