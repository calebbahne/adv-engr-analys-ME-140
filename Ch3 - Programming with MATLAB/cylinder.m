function cylinder(r,L,plot_title)
%cylinder: problem 3.11 plots volume of liquid in a horiz cylinder
% inputs
%   r = radius
%   L = length
%   plot_title

h = linspace(0,2*r);
V = (r.^2.*acos((r-h)./r)-(r-h).*sqrt(2*r.*h-h.^2)).*L;
plot(h,V);
title(plot_title);
end