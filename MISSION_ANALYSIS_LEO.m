% PRELIMINARY MISSION ANALYSIS OF A SATELLITE IN LOW EARTH ORBIT (LEO).
% Requested inputs: desired satellite altitude and the latitude/longitude 
% of the observation point on Earth. The input satellite will follow 
% a sun-synchronous orbit (SSO).
% If a specific revisit time is needed (e.g., passing over a target along 
% the GT once every 2 days), the corresponding SS-O altitude can be selected.

clear all

% Input definition
prompt = {'Enter satellite altitude h (in km):','Enter observation point latitude:','Enter observation point longitude:'};
dlgtitle = 'Input';
dims = [1 50];
definput = {'621.86','-75','-123'};
opts.Interpreter = 'tex';
answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
h = str2num(answer{1});
lat_TARGET = str2num(answer{2});
lon_TARGET = str2num(answer{3});

% Input validation
if h<500 || h>1000
    error('h must be between 500km and 1000km (LEO orbit requires h < 1000km, while below 500km atmospheric DRAG effect is excessive)')
end
if lat_TARGET<-80 || lat_TARGET>80 || lon_TARGET<-180 || lon_TARGET>180
    error('Target latitude must be between -80° and 80°, and longitude between -180° and 180°')
end

% Initial parameters definition
mu=398600.44;
J2=0.00108263;
R_earth=6378.14;
a=R_earth+h;
a=round(a*100)/100;
e=0;
RAAN=0;
omega=0;
ni=0;
nn=sqrt(mu/(a)^3);
rev_rate=(2*pi)/(86400*365.242199);
i=acos(rev_rate/((-3/2)*J2*((R_earth/a)^2)*nn));
i=round(i*100)/100;
i_deg=rad2deg(i);

% Conversion from initial orbital parameters to initial position and velocity
% + integration of the restricted two-body problem equation
[r0,v0] = COE2RV(a,e,i,RAAN,omega,ni);
tspan=0:10:24*60*60;
y=[r0;v0];
options=odeset('RelTol',10^(-12),'AbsTol', 10^(-12));
[t,y]=ode113(@twobodyproblem,tspan,y,options);
r=y(:,1:3);
v=y(:,4:6);

% Orbit visualization
t_REV=2*pi*sqrt(a^3/mu);
rev_instants = (t_REV*8641/(60*60*24))-7;
figure(1)
plot3(r(:,1),r(:,2),r(:,3),'r.');
hold on
[x,y,z]=ellipsoid(0,0,0,6378.14,6378.14,6378.14,128);
surf(x,y,-z);
axis equal
grid on
title('Satellite orbit')

% Ground track visualization
N=size(r,1);
G0=0;
ne=360/86164;
for j=1:N
    lat(j)=asind(r(j,3)/norm(r(j,:)));
    lon_gamma(j)=atan2d(r(j,2)/(norm(r(j,:))*cosd(lat(j))),r(j,1)/(norm(r(j,:))*cosd(lat(j))));
    G(j)=G0+ne*t(j);
    lon(j)=lon_gamma(j)-G(j);
    n_revs=abs(round(lon(j)/360));
    while lon(j)<-180
        lon(j)=lon(j)+360;
    end
    while lon(j)>180
        lon(j)=lon(j)-360;
    end
end
hold off
figure(2)
plot(lon,lat,'r.')
xlim([-180 180]);
ylim([-90 90]);
xlabel('Longitude');
ylabel('Latitude');
title('Ground Track');
hold on
I=imread('Earth.jpg');
H=image(xlim,-ylim,I);
uistack(H,'bottom');
legend('Ground Track');

% Visibility
angspeed_earth = 2*pi/86164;
x_ECI = R_earth.*cosd(lat_TARGET).*cosd(lon_TARGET);
y_ECI = R_earth.*cosd(lat_TARGET).*sind(lon_TARGET);
radius_target = norm([x_ECI y_ECI 0]);
x_TARGET = radius_target.*cos(angspeed_earth.*(tspan')+deg2rad(lon_TARGET)); 
y_TARGET = radius_target.*sin(angspeed_earth.*(tspan')+deg2rad(lon_TARGET));
z_TARGET = ones(size(x_TARGET,1),1).*R_earth.*sind(lat_TARGET);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot3(r(:,1),r(:,2),r(:,3),'r.');
hold on
plot3(x_TARGET, y_TARGET, z_TARGET,'b.');
[x,y,z]=ellipsoid(0,0,0,6378.14,6378.14,6378.14,128);
surf(x,y,-z);
axis equal
grid on
title('Satellite orbit')
legend('Satellite orbit','Ground station')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

target_vector=[x_TARGET, y_TARGET, z_TARGET];
target_unit_vector = target_vector./norm(target_vector);
rho=r-target_vector;
rho_unit_vector = rho./norm(rho);
rho_max = R_earth*(sqrt(((h+R_earth)/R_earth)^2-(cosd(5))^2)-sind(5));

% Count visibility instances
counter=1;
visibility_instants=[];
while(counter<=length(tspan))
    if norm(rho(counter,:))<=rho_max
        visibility_instants(counter)=1;
    else
        visibility_instants(counter)=0;
    end
    counter=counter+1;
end
sum_instants=sum(visibility_instants);
visibility_time=(sum_instants/length(tspan))*24*60;

% Visibility check
counter=1;
EL=[];
visibility_instants_check=[];
while(counter<=length(tspan))
    EL(counter)=asin(((h*(h+2*R_earth))-(norm(rho(counter,:)))^2)/(2*R_earth*norm(rho(counter,:))));
    if EL(counter)>=deg2rad(5)
        visibility_instants_check(counter)=1;
    else
        visibility_instants_check(counter)=0;
    end
    counter=counter+1;
end
sum_instants=sum(visibility_instants_check);
visibility_time_check=(sum_instants/length(tspan))*24*60;

% Conversion to minutes and seconds
minutes_val = floor(visibility_time);
seconds_val = round((visibility_time-floor(visibility_time))*60);

% Print Output
fprintf('Visibility time is %3d minutes and %2d seconds within a 24-hour interval.', minutes_val, seconds_val)
msgbox("You can read the satellite visibility time in the Command Window.","Analysis Complete")