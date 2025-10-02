function Gdir=G_dir_function(ng,xg,wg,gL,hL, muprime,phiprime)
% this routine calculates the G_FUNCTION given a direction of
% photon travel OMEGA^PRIME and the leaf normal pdf

%------Variables ------
%   xg(ng), wg(ng), gL(ng), hL(ng)
%   muprime, phiprime
% 	mu_tp, sin_tp
% 	upperlimit_tL, lowerlimit_tL
%   conv1_tL, conv2_tL, sum_tL
% 	upperlimit_pL, lowerlimit_pL
%   conv1_pL, conv2_pL, sum_pL
%  	neword_tL, mu_tL, sin_tL
% 	neword_pL, dotproduct, Gdir
% end declarations

mu_tp  = muprime;
sin_tp = sqrt(1.0 - muprime*muprime);

% define limits of integration and the convertion factors for integration
% over thetaL (note the tL suffix!)
upperlimit_tL = pi/2.0;
lowerlimit_tL = 0.0;
conv1_tL = (upperlimit_tL-lowerlimit_tL)/2.0;
conv2_tL = (upperlimit_tL+lowerlimit_tL)/2.0;

% define limits of integration and the convertion factors for integration
% over phiL (note the pL suffix!)
upperlimit_pL = 2.0*pi;
lowerlimit_pL = 0.0;
conv1_pL = (upperlimit_pL-lowerlimit_pL)/2.0;
conv2_pL = (upperlimit_pL+lowerlimit_pL)/2.0;

% integral over theta_L, the outer layer loop
neword_tL = conv1_tL*xg+conv2_tL;
mu_tL     = cos(neword_tL'); mu_tL = repmat(mu_tL,[ng 1]);
sin_tL    = sin(neword_tL'); sin_tL = repmat(sin_tL,[ng 1]);

% integral over phi_L, the inner layer loop
neword_pL = conv1_pL*xg + conv2_pL;
dotproduct = abs (mu_tL*mu_tp +sin_tL*sin_tp.*cos(neword_pL-phiprime));
sum_pL = sum(wg.*hL./(2.0*pi).*dotproduct,1);
%end of the inner layer loop

% finish the phi_L integral
sum_pL = sum_pL*conv1_pL;
sum_tL = sum(wg.*gL.*sum_pL');
% end of the outer layer loop

% finish the theta_L integral
sum_tL = sum_tL*conv1_tL;
Gdir  = sum_tL;
end

