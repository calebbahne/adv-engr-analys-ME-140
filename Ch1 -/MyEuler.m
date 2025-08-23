function [t,y] = MyEuler(f,a,b,N,alp)
% Does Euler's Method to calculate y(t) for an IVP. 
% from f(t, y, y', y'', ..., y^n) = 0
% f =@(t,y) = y' 
% N is the number of sections you're dividing a->b into
% alp = y(t0)

% Example input for y' = t + y
%    f=@(x,y) x + y; a =1; b=2; N =10; alp=1;
%    [t,y]=MyEuler(f,a,b,N,alp);
%           g = @(t) -t-1+3/exp(1)*exp(t)
%    plot(t,y);  hold on; plot(t,g(t));

h=(b-a)/N;
t=a:h:b; %either a range or x=linspace(a,b,N+1)
y=zeros(1,N+1);
y(1)=alp;
for i=2:N+1
    y(i)=y(i-1)+h*f(t(i-1),y(i-1));
end