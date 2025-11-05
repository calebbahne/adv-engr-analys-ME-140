Tinf = 250;
Tmi = 1600;
Tmo = 1400;
mdot = 3.647;
cp = 1.087e3;
Di = 1; L = 1; kins = 0.125;
hbar = 95.02;

(Tinf-Tmo)/(Tinf-Tmi)

f = @(x) exp(-1/(mdot*cp*( log((Di/2+x)/(Di/2))/(2*pi*L*kins) + 1/(hbar*pi*Di*L) )))  - ((Tinf-Tmo)/(Tinf-Tmi));
% help fzero
fzero(f, 0.005)
fplot(f, [-0.02 0.01])

A = -1/log(0.8519)/mdot/cp

t = exp(0.7854*(A-1/(hbar*pi*Di*L)))*0.5-0.5