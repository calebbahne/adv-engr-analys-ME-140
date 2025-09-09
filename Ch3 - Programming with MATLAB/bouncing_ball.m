%Problem 3.21 & 2.
clc,clf,clear
maxit=1000;
g=9.81; theta0=50*pi/180; v0=5; CR=0.83;
j=1;t(j)=0;x=0;y=0;
xx=x;yy=y;
plot(x,y,'o','MarkerFaceColor','b','MarkerSize',8)
xmax=8; axis([0 xmax 0 0.8])
M(1)=getframe;
dt=1/128;
j=1; xxx=0; iter=0;
while(1)
tt=0;
timpact=2*v0*sin(theta0)/g;
ximpact=v0*cos(theta0)*timpact;
while(1)
j=j+1;
h=dt;
if tt+h>timpact,h=timpact-tt;end
t(j)=t(j-1)+h;
tt=tt+h;
x=xxx+v0*cos(theta0)*tt;
y=v0*sin(theta0)*tt-0.5*g*tt^2;
xx=[xx x];yy=[yy y];
plot(xx,yy,':',x,y,'o','MarkerFaceColor','b','MarkerSize',8)
axis([0 xmax 0 0.8])
M(j)=getframe;
iter=iter+1;
if tt>=timpact, break, end
end
v0=CR*v0;
xxx=x;
if x>=xmax|iter>=maxit,break,end
end
pause
clf
axis([0 xmax 0 0.8])
movie(M,1,36)
