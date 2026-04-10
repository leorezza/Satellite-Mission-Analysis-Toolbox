% PRELIMINARY MISSION ANALYSIS OF SATELLITES IN LOW EARTH ORBIT (LEO):
% CONSTELLATION of 4 satellites in SS-O orbits (whose initial parameters 
% differ only by the RAAN value).
%
% This script, using the function MISSION_ANALYSIS_FUNCTION, extracts 
% vectors with components equal to zero or one depending on whether the 
% satellite is visible from an observation point at that specific time instant. 
% To define the satellite orbit and the observation point, you must provide 
% the desired SS-O altitude and the coordinates (latitude and longitude) 
% of the observation point – additionally, you must enter the RAAN of the 
% orbit (in radians!).
% The script performs a cross-sum of the visibility instants of the 
% entered satellites, providing an estimate of the total visibility time 
% over a day for the entire constellation.
% If you wish to fly over a specific point with a certain revisit time 
% (e.g., passing over a target once every 2 days), you can select the 
% corresponding SS-O altitude (e.g., 2D29R).

% Input definitions
prompt = {'Enter satellite altitude h (in km):','Enter observation point latitude:','Enter observation point longitude:','Enter RAAN satellite 1:','Enter RAAN satellite 2:','Enter RAAN satellite 3:','Enter RAAN satellite 4:','M1:','M2:','M3:','M4'};
dlgtitle = 'Input';
dims = [1 50];
definput = {'621.86','-75','-123','0','90*pi/180','180*pi/180','270*pi/180','0','90*pi/180','pi','270*pi/180'};
opts.Interpreter = 'tex';
answer = inputdlg(prompt,dlgtitle,dims,definput,opts);
h = str2num(answer{1});
lat_obs = str2num(answer{2});
lon_obs = str2num(answer{3});
RAAN1 = str2num(answer{4});
RAAN2 = str2num(answer{5});
RAAN3 = str2num(answer{6});
RAAN4 = str2num(answer{7});
M1 = str2num(answer{8});
M2 = str2num(answer{9});
M3 = str2num(answer{10});
M4 = str2num(answer{11});

% Visibility instants of the constellation
[I_satellite1,r1,t] = MISSION_ANALYSIS_FUNCTION(h,lat_obs,lon_obs,RAAN1,M1);
[I_satellite2,r2] = MISSION_ANALYSIS_FUNCTION(h,lat_obs,lon_obs,RAAN2,M2);
[I_satellite3,r3] = MISSION_ANALYSIS_FUNCTION(h,lat_obs,lon_obs,RAAN3,M3);
[I_satellite4,r4] = MISSION_ANALYSIS_FUNCTION(h,lat_obs,lon_obs,RAAN4,M4);

% Visualization of satellite ground tracks
N=size(r1,1);
G0=0;
ne=360/86164;

% Ground track calculation for Satellite 1
for j=1:N
    lat1(j)=asind(r1(j,3)/norm(r1(j,:)));
    lon_gamma1(j)=atan2d(r1(j,2)/(norm(r1(j,:))*cosd(lat1(j))),r1(j,1)/(norm(r1(j,:))*cosd(lat1(j))));
    G(j)=G0+ne*t(j);
    lon1(j)=lon_gamma1(j)-G(j);
    n_revs=abs(round(lon1(j)/360));
    while lon1(j)<-180
        lon1(j)=lon1(j)+360;
    end
    while lon1(j)>180
        lon1(j)=lon1(j)-360;
    end
end

% Ground track calculation for Satellite 2
for j=1:N
    lat2(j)=asind(r2(j,3)/norm(r2(j,:)));
    lon_gamma2(j)=atan2d(r2(j,2)/(norm(r2(j,:))*cosd(lat2(j))),r2(j,1)/(norm(r2(j,:))*cosd(lat2(j))));
    G(j)=G0+ne*t(j);
    lon2(j)=lon_gamma2(j)-G(j);
    while lon2(j)<-180
        lon2(j)=lon2(j)+360;
    end
    while lon2(j)>180
        lon2(j)=lon2(j)-360;
    end
end

% Ground track calculation for Satellite 3
for j=1:N
    lat3(j)=asind(r3(j,3)/norm(r3(j,:)));
    lon_gamma3(j)=atan2d(r3(j,2)/(norm(r3(j,:))*cosd(lat3(j))),r3(j,1)/(norm(r3(j,:))*cosd(lat3(j))));
    G(j)=G0+ne*t(j);
    lon3(j)=lon_gamma3(j)-G(j);
    while lon3(j)<-180
        lon3(j)=lon3(j)+360;
    end
    while lon3(j)>180
        lon3(j)=lon3(j)-360;
    end
end

% Ground track calculation for Satellite 4
for j=1:N
    lat4(j)=asind(r4(j,3)/norm(r4(j,:)));
    lon_gamma4(j)=atan2d(r4(j,2)/(norm(r4(j,:))*cosd(lat4(j))),r4(j,1)/(norm(r4(j,:))*cosd(lat4(j))));
    G(j)=G0+ne*t(j);
    lon4(j)=lon_gamma4(j)-G(j);
    while lon4(j)<-180
        lon4(j)=lon4(j)+360;
    end
    while lon4(j)>180
        lon4(j)=lon4(j)-360;
    end
end

hold off
figure(1)
plot(lon1,lat1,'r.')
hold on
plot(lon2,lat2,'b.')
hold on
plot(lon3,lat3,'g.')
hold on
plot(lon4,lat4,'y.')
xlim([-180 180]);
ylim([-90 90]);
xlabel('Longitude');
ylabel('Latitude');
title('Ground Tracks');
hold on
I_map=imread('Earth.jpg');
H_image=image(xlim,-ylim,I_map);
uistack(H_image,'bottom');

% Total Visibility
I_total = I_satellite1+I_satellite2+I_satellite3+I_satellite4;
for i=1:8641
    if I_total(i)>1
        I_total(i)=1;
    end
end

% Time Conversion
seconds_total = sum(I_total)*86126/8641;
minutes_total = seconds_total/60;
hours_total = minutes_total/60;

hou = floor(hours_total);
minutes_rem = (hours_total-hou)*60;
min_val = floor(minutes_rem);
seconds_rem = (minutes_rem-floor(minutes_rem))*60;
sec_val = round(seconds_rem);

% Print results
fprintf('The constellation visibility time is %2d hours, %2d minutes and %2d seconds - within a 24-hour interval.',hou,min_val,sec_val)
msgbox("You can read the constellation visibility time in the Command Window.","Analysis Complete")