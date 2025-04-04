function az = simpleAzimuth(lon1, lat1, lon2, lat2)


% calculate the simple backazimuth from station

R = 6371000; % m

phi1 = deg2rad(lat1);
phi2 = deg2rad(lat2);
lambda1 = deg2rad(lon1);
lambda2 = deg2rad(lon2);
dphi = deg2rad(lat2 - lat1);
dlambda = deg2rad(lon2 - lon1);

% a = (sin(dphi/2))^2 +...
%     cos(phi1) .* cos(phi2) .*...
%     (sin(dlambda/2)).^2;

%midpoint
% c = 2 .* atan2(sqrt(a), sqrt(1-a));

% distance between points
% d = R .* c;

% bearing:
y = sin(lambda2 - lambda1) .* cos(phi2);
x = cos(phi1).*sin(phi2) - sin(phi1).*cos(phi2).*cos(lambda2 - lambda1);
az = atan2(y, x);

az = rad2deg(az);
% baz = a;