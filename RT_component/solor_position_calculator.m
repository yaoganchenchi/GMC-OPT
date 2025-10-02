function [sunzenith,sunazimuth]=solor_position_calculator(iyear,idoy,lon,lat,localtime)
% https://gml.noaa.gov/grad/solcalc/solareqns.PDF

%% test
% iyear = 2021;
% idoy = 243;
% lon = -71.05; % Boston
% lat = 42.35;
% 
% lon = -122.416777; % San Fran
% lat = 37.77777;
% localtime = 9.5;

rad2deg = 180.0/pi;
deg2rad = pi/180.0;

%% fractional year gamma
if ~isLeapYear(iyear)
    gamma = 2*pi/365*(idoy-1+ (localtime-12)/24);
else
    gamma = 2*pi/366*(idoy-1+ (localtime-12)/24);
end

%% equation of time and solar declination angle
eqtime = 229.18*(0.000075 + 0.001868*cos(gamma) - 0.032077*sin(gamma) - 0.014615*cos(2*gamma)-0.040849*sin(2*gamma)); % in minutes
decl = 0.006918 - 0.399912*cos(gamma) + 0.070257*sin(gamma) - 0.006758*cos(2*gamma) + 0.000907*sin(2*gamma) - 0.002697*cos(3*gamma) + 0.00148*sin (3*gamma);

%% time offset
timezone = round(lon/15); % in hour
time_offset = eqtime + 4*lon - 60*timezone;

%% true solar time
hr = floor(localtime);
mn = floor(60*(localtime - floor(localtime)));
sc = 60*(60*(localtime - floor(localtime))-floor(60*(localtime - floor(localtime))));
tst = hr*60 + mn + sc/60 + time_offset;

%% solar hour angle
ha = (tst/4)-180;

%% calculate solar zenith angle
cos_sunzenith = sin(decl)*sin(lat*deg2rad) + cos(decl)*cos(lat*deg2rad)*cos(ha*deg2rad);
sunzenith = acos(cos_sunzenith)*rad2deg;

%% calculate solar azimuth angle
num = cos(decl)*cos(ha*deg2rad); % numerator 
den = (cos(decl)*cos(lat*deg2rad)*cos(ha*deg2rad)-sin(decl)*cos(lat*deg2rad)); % demominator 
tan_sunazimuth = num / den; 
if num>0 && den>0
    sunazimuth = 180-atan(tan_sunazimuth)*rad2deg;
elseif num>0 && den<0
    sunazimuth = 0-atan(tan_sunazimuth)*rad2deg;
elseif num<0 && den>0
    sunazimuth = 180-atan(tan_sunazimuth)*rad2deg;
elseif num<0 && den<0
    sunazimuth = 360-atan(tan_sunazimuth)*rad2deg;
end

end

