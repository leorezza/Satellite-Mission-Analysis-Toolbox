% Transformation from orbital parameters to initial position and velocity

function [r,v] = COE2RV(a,e,i,raan,omega,ni)
mu=4902.7779; % Gravitational parameter of the Moon
p=a*(1-e^2);
h=sqrt(mu*p);
e_v=[cos(omega)*cos(raan)-sin(omega)*cos(i)*sin(raan);cos(omega)*sin(raan)+sin(omega)*cos(i)*cos(raan);sin(omega)*sin(i)];
h_v=[sin(i)*sin(raan); -sin(i)*cos(raan); cos(i)];
p_v=cross(h_v,e_v);

r=p/(1+e*cos(ni))*(cos(ni)*e_v+sin(ni)*p_v);
v=mu/h*(e*p_v+(-sin(ni)*e_v+cos(ni)*p_v));
