% Matlab numerical solution using Euler's method from Math 121
g = 9.81; % kg*m/s^2
m = 68.1; c_d = 0.25; % kg
a = 0; b = 12; alp = 0;
f = @(t,y) g - c_d/m*y^2;

% Step size 1
h=1; % step size
t=a:h:b; %either a range or x=linspace(a,b,N+1)
N=12;
y1=zeros(1,N+1);
y1(1)=alp;
for i=2:N+1
    y1(i)=y1(i-1)+h*f(t(i-1),y1(i-1));
end
y1(N+1);

% Step size 0.2
h=0.2; % step size
t2 = a:h:b;
N=12*5;
y0_2=zeros(1,N+1);
y0_2(1)=alp;
for i=2:N+1
    y0_2(i)=y0_2(i-1)+h*f(t2(i-1),y0_2(i-1));
end
y0_2(N+1);

% Plot results
plot(t,y1, 'r');
hold on;
plot(t2,y0_2, 'g');
v_true=zeros(1,N+1);
v_true = sqrt(g*m/c_d)*tanh(sqrt(g*c_d/m)*t);
plot(t,v_true, 'b');
title('Bungee Jumper Velocity');
xlabel('Time (sec)');
ylabel('Velocity (m/s)');
legend('Num 1','Num 0.2','Analytical','Location','southeast');
hold off;