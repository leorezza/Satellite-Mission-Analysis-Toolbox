% PRELIMINARY MISSION ANALYSIS OF SATELLITES IN LOW EARTH ORBIT (LEO).
% The function requires as input the desired satellite altitude, the 
% latitude and longitude of the target to be overflown, and the RAAN 
% of the desired orbit.
% If a specific revisit time is needed (e.g., passing over the target 
% once every 2 days), the corresponding SS-O altitude can be selected, 
% e.g., 2D29R, by first calculating the nodal orbital period P = 86400/R, 
% then the semi-major axis a = (mu*(P/2*pi)^2)^1/3 and finally h = a - R_earth. 
% This provides the altitude to be used as input.

function [visibility_instants, r, t] = MISSION_ANALYSIS_FUNCTION(h, lat_TARGET, lon_TARGET, RAAN, M)
% Input validation
if h < 500 || h > 1000
    error('h must be between 500km and 1000km (LEO orbit requires h < 1000km, while below 500km atmospheric DRAG effect is excessive)')
end
if lat_TARGET < -80 || lat_TARGET > 80 || lon_TARGET < -180 || lon_TARGET > 180
    error('Target latitude must be between -80° and 80°, while longitude must be between -180° and 180°')
end

% Initial parameters definition
mu = 398600.44;
J2 = 0.00108263;
R_earth = 6378.14;
a = R_earth + h;
a = round(a * 100) / 100;
e = 0;
omega = 0;
ni = 0;
nn = sqrt(mu / (a)^3);
rev_rate = (2 * pi) / (86400 * 365.242199);
i = acos(rev_rate / ((-3/2) * J2 * ((R_earth / a)^2) * nn));
i = round(i * 100) / 100;
i_deg = rad2deg(i);

% Conversion from initial orbital parameters to initial position and velocity
% + integration of the restricted two-body problem equation
[r0, v0] = COE2RV(a, e, i, RAAN, omega, M);
tspan = 0:10:24*60*60;
y = [r0; v0];
options = odeset('RelTol', 10^(-12), 'AbsTol', 10^(-12));
[t, y] = ode113(@twobodyproblem, tspan, y, options);
r = y(:, 1:3);
v = y(:, 4:6);

% GROUND TRACK VISUALIZATION
N = size(r, 1);
G0 = 0;
ne = 360 / 86164;
for j = 1:N
    lat_val(j) = asind(r(j, 3) / norm(r(j, :)));
    lon_gamma(j) = atan2d(r(j, 2) / (norm(r(j, :)) * cosd(lat_val(j))), r(j, 1) / (norm(r(j, :)) * cosd(lat_val(j))));
    G(j) = G0 + ne * t(j);
    lon_val(j) = lon_gamma(j) - G(j);
    n_revs = abs(round(lon_val(j) / 360));
    while lon_val(j) < -180
        lon_val(j) = lon_val(j) + 360;
    end
    while lon_val(j) > 180
        lon_val(j) = lon_val(j) - 360;
    end
end

hold off
figure(1)
plot(lon_val, lat_val, 'r.')
xlim([-180 180]);
ylim([-90 90]);
xlabel('Longitude');
ylabel('Latitude');
title('Ground Track');
hold on
I = imread('Earth.jpg');
H = image(xlim, -ylim, I);
uistack(H, 'bottom');

%==========================================================================
% Visibility
angspeed_earth = 2 * pi / 86164;
x_ECI = R_earth .* cosd(lat_TARGET) .* cosd(lon_TARGET);
y_ECI = R_earth .* cosd(lat_TARGET) .* sind(lon_TARGET);
radius_target = norm([x_ECI y_ECI 0]);

x_TARGET = radius_target .* cos(angspeed_earth .* (tspan') + deg2rad(lon_TARGET)); 
y_TARGET = radius_target .* sin(angspeed_earth .* (tspan') + deg2rad(lon_TARGET));
z_TARGET = ones(size(x_TARGET, 1), 1) .* R_earth .* sind(lat_TARGET);

target_vector = [x_TARGET, y_TARGET, z_TARGET];
% Note: The unit vector calculation remains but is not actively used in the loop below
target_unit_vector = target_vector ./ norm(target_vector); 

rho = r - target_vector;
% Note: The rho unit vector calculation remains as per original logic
rho_unit_vector = rho ./ norm(rho);

rho_max = R_earth * (sqrt(((h + R_earth) / R_earth)^2 - (cosd(5))^2) - sind(5));

% Count visibility instances
counter = 1;
visibility_instants = [];
while(counter <= length(tspan))
    if norm(rho(counter, :)) <= rho_max
        visibility_instants(counter) = 1;
    else
        visibility_instants(counter) = 0;
    end
    counter = counter + 1;
end

sum_instants = sum(visibility_instants);
visibility_time = (sum_instants / length(tspan)) * 24 * 60;