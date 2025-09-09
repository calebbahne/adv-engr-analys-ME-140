function phasor(r, nt, nm)
% function to show the orbit of a phasor
% r = radius
% nt = number of increments for theta
% nm = number of movies
clc;clf
dtheta=2*pi/nt;
th=0;
fac=1.2;
xx=r;yy=0;
for i=1:nt+1
x=r*cos(th);y=r*sin(th);
xx=[xx x];yy=[yy y];
plot([0 x],[0 y],xx,yy,':',...
x,y,'o','MarkerFaceColor','b','MarkerSize',8)
axis([-fac*r fac*r -fac*r fac*r]);
axis square
M(i)=getframe;
th=th+dtheta;
end
pause
clf
axis([-fac*r fac*r -fac*r fac*r]);
axis square
movie(M,1,36)
