function [gL,hL]=leaf_normal_pdf(xg,wg,gL,hL,lnd)
% this routine evaluates PLANOPHILE leaf normal inclination pdf (gL)
% and UNIFORM leaf normal azimuthal pdf (hL)

%------Variables ------
%  xg(ng), wg(ng), gL(ng), hL(ng)
%  upperlimit, lowerlimit, conv1, conv2, sum
%  neword

hL(:) = 1.0;
% obtain the planophile gL
%   gL(thetaL) = (2/pi) (1+cos(2thetaL))
% and at the same time check if it satisfies the required condition
% of normalization, that is,
%   [ int_0^(pi/2) dthetaL gL(thetaL) = 1.0 ]

%  step 1: define limits of integration and the convertion factors
upperlimit = pi/2.0;
lowerlimit = 0.0;
conv1 = (upperlimit-lowerlimit)/2.0;
conv2 = (upperlimit+lowerlimit)/2.0;

%  step 2: do the integral by making sure the ordinates run between
%          the upperlimit and the lowerlimit
neword = conv1*xg + conv2;
if strcmp(lnd,'Planophile') % mostly horizontal
    gL = 2.0/pi*(1+cos(2.0*neword));
elseif strcmp(lnd,'Erectophile') % mostly vertical
    gL = 2.0/pi*(1-cos(2.0*neword));
elseif  strcmp(lnd,'Plagiophile') % mostly 45 deg
    gL = 2.0/pi*(1-cos(4.0*neword));
elseif  strcmp(lnd,'Extremophile') % half horizontal, half vertical
    gL = 2.0/pi*(1+cos(4.0*neword));
elseif  strcmp(lnd,'Uniform') % Any angle is possible
    gL = 2.0/pi;
elseif strcmp(lnd,'Spherical')
    gL = sin(neword);
end

%  step 3: make sure not to forget to apply the conversion factor 1 again
%{
sum0 = sum(gL.*wg);
sum0 = sum0*conv1;
fprintf('LNO  check (=1.0?): %.8f\n',sum0);
%}

end

