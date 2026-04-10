clear all; clc;

% Initial parameters
h=100; % orbit altitude of LIGMA A and A'
h2=80; % orbit altitude of LIGMA B and C
mu=4902.7779; % Gravitational parameter of the Moon
J2=0.0002027; % J2 factor
R=1737.1; % Moon's gravitational parameter
a=R+h; % Semiaxis major of LIGMA A and A'
a2=R+h2; % Semiaxis majorof LIGMA B and C
e=0; % eccentricity
RAAN=0; % RAAN of LIGMA A, A' an B
RAAN2=1; % RAAN of LIGMA A1
omega=0; % Argument of periapsis
ni=0; % True anomaly of LIGMA A and B
ni2=6.25*pi/180; % True anomaly of LIGMA A'
ni3=90*pi/180; % True anomaly of LIGMA C
i=89*pi/180; % Inclination of LIGMA A and A'
i2=87.4*pi/180; % Inclination of LIGMA C and C
% T=2*pi*sqrt((a)^3/mu); % Orbital period

% Conversion from initial orbital parameters to initial position and velocity
% and integration of the restricted two-body problem equation
% LIGMA A
[r0,v0] = COE2RV(a,e,i,RAAN,omega,ni);
tspan=0:0.5:2*60*60;
y=[r0;v0];
settings=odeset('RelTol',10^(-12),'AbsTol', 10^(-12));
[t,y]=ode113(@twobodyproblem,tspan,y,settings);
r=y(:,1:3);
v=y(:,4:6);
% LIGMA A'
[r0,v0] = COE2RV(a,e,i,RAAN,omega,ni2);
tspan=0:0.5:2*60*60;
y=[r0;v0];
settings=odeset('RelTol',10^(-12),'AbsTol', 10^(-12));
[t,y]=ode113(@twobodyproblem,tspan,y,settings);
r2=y(:,1:3);
v3=y(:,4:6);
% LIGMA B
[r0,v0] = COE2RV(a2,e,i2,RAAN,omega,ni);
tspan=0:0.5:2*60*60;
y=[r0;v0];
settings=odeset('RelTol',10^(-12),'AbsTol', 10^(-12));
[t,y]=ode113(@twobodyproblem,tspan,y,settings);
r3=y(:,1:3);
v3=y(:,4:6);
% LIGMA C
[r0,v0] = COE2RV(a2,e,i2,RAAN2,omega,ni3);
tspan=0:0.5:2*60*60;
y=[r0;v0];
settings=odeset('RelTol',10^(-12),'AbsTol', 10^(-12));
[t,y]=ode113(@twobodyproblem,tspan,y,settings);
r4=y(:,1:3);
v4=y(:,4:6);

% Orbits 3D visualiation
t_REV=2*pi*sqrt(a^3/mu);
ist_REV=(t_REV/(29.5*60*60*24));
figure(1)
plot3(r(:,1),r(:,2),r(:,3),'k.');
hold on
plot3(r2(:,1),r2(:,2),r2(:,3),'r.');
plot3(r3(:,1),r3(:,2),r3(:,3),'g.');
plot3(r4(:,1),r4(:,2),r4(:,3),'b.');
[x,y,z]=ellipsoid(0,0,0,1738.1,1738.1,1736.0,128);
h=surf(x,y,-z);
textureImage = imread('moonmap.jpg');
set(h, 'FaceColor', 'texturemap', 'EdgeColor', 'none', 'CData', textureImage);
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on
title('Satellites orbits')
legend('LIGMA A','LIGMA A''','LIGMA B','LIGMA C');

% Ground tracks visualiation
N=size(r,1);
G0=0;
ne=360/(27.3*24*60*60);
for j=1:N
    lat(j)=asind(r(j,3)/norm(r(j,:)));
    lon_gamma(j)=atan2d(r(j,2)/(norm(r(j,:))*cosd(lat(j))),r(j,1)/(norm(r(j,:))*cosd(lat(j))));
    G(j)=G0+ne*t(j);
    lon(j)=lon_gamma(j)-G(j);
    n_giri=abs(round(lon(j)/360));
    while lon(j)<-180
        lon(j)=lon(j)+360;
    end
    while lon(j)>180
        lon(j)=lon(j)-360;
    end
end
for j=1:N
    lat2(j)=asind(r2(j,3)/norm(r2(j,:)));
    lon_gamma2(j)=atan2d(r2(j,2)/(norm(r2(j,:))*cosd(lat2(j))),r2(j,1)/(norm(r2(j,:))*cosd(lat2(j))));
    G(j)=G0+ne*t(j);
    lon2(j)=lon_gamma2(j)-G(j);
    n_giri=abs(round(lon2(j)/360));
    while lon2(j)<-180
        lon2(j)=lon2(j)+360;
    end
    while lon2(j)>180
        lon2(j)=lon2(j)-360;
    end
end
for j=1:N
    lat3(j)=asind(r3(j,3)/norm(r3(j,:)));
    lon_gamma3(j)=atan2d(r3(j,2)/(norm(r3(j,:))*cosd(lat3(j))),r3(j,1)/(norm(r3(j,:))*cosd(lat3(j))));
    G(j)=G0+ne*t(j);
    lon3(j)=lon_gamma3(j)-G(j);
    n_giri=abs(round(lon3(j)/360));
    while lon3(j)<-180
        lon3(j)=lon3(j)+360;
    end
    while lon3(j)>180
        lon3(j)=lon3(j)-360;
    end
end
for j=1:N
    lat4(j)=asind(r4(j,3)/norm(r4(j,:)));
    lon_gamma4(j)=atan2d(r4(j,2)/(norm(r4(j,:))*cosd(lat4(j))),r4(j,1)/(norm(r4(j,:))*cosd(lat4(j))));
    G(j)=G0+ne*t(j);
    lon4(j)=lon_gamma4(j)-G(j);
    n_giri=abs(round(lon4(j)/360));
    while lon4(j)<-180
        lon4(j)=lon4(j)+360;
    end
    while lon4(j)>180
        lon4(j)=lon4(j)-360;
    end
end
hold off
figure(2)
plot(lon,lat,'k.')
hold on
plot(lon2,lat2,'r.')
hold on
plot(lon3,lat3,'g.')
hold on
plot(lon4,lat4,'b.')
xlim([-180 180]);
ylim([-90 90]);
xlabel('Longitude');
ylabel('Latitude');
title('Ground Tracks');
hold on
I=imread('Moon.jpg');
H=image(xlim,-ylim,I);
uistack(H,'bottom');
legend('LIGMA A','LIGMA A''','LIGMA B','LIGMA C');