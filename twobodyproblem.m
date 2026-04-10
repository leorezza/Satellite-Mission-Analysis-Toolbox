% Restricted two-body problem equation

function [Ydot]=twobodyproblem(t,Y)

mu=4902.7779; % Gravitational parameter of the Moon
x=Y(1);
y=Y(2);
z=Y(3);
vx=Y(4);
vy=Y(5);
vz=Y(6);

r=norm(Y(1:3));

xdot=vx;
ydot=vy;
zdot=vz;

vxdot=-mu/r^3*x;
vydot=-mu/r^3*y;
vzdot=-mu/r^3*z;

Ydot=[xdot ydot zdot vxdot vydot vzdot]';
