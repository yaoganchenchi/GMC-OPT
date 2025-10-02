function Gamma_d=GAMMA_d_function (ng,xg,wg,gL,hL,rho_Ld,tau_Ld, muprime,phiprime,mu,phi)
% this routine calculates the G_FUNCTION given a direction of
% photon travel OMEGA^PRIME and the leaf normal pdf

%------Variables ------
%   xg(ng), wg(ng), gL(ng), hL(ng), rho_Ld, tau_Ld
%   muprime, phiprime, mu, phi
% 	mu_t, mu_tp, sin_t, sin_tp
% 	upperlimit_tL, lowerlimit_tL
%   conv1_tL, conv2_tL, sum_tL
% 	upperlimit_pL, lowerlimit_pL
%   conv1_pL, conv2_pL, sum_pL
%  	neword_tL, mu_tL, sin_tL
% 	neword_pL, dotproduct1, dotproduct2, Gamma_d

mu_t   = mu;
mu_tp  = muprime;
sin_t  = sqrt(1.0 - mu*mu);
sin_tp = sqrt(1.0 - muprime*muprime);

% define limits of integration and the convertion factors for integration
% over thetaL (note the tL suffix!)

upperlimit_tL = pi/2.0;   % why pi/2? Leaf normal is hemispherical, not spherical
lowerlimit_tL = 0.0;
conv1_tL = (upperlimit_tL-lowerlimit_tL)/2.0;
conv2_tL = (upperlimit_tL+lowerlimit_tL)/2.0;

% define limits of integration and the convertion factors for integration
% over phiL (note the pL suffix!)

upperlimit_pL = 2.0*pi;
lowerlimit_pL = 0.0;
conv1_pL = (upperlimit_pL-lowerlimit_pL)/2.0;
conv2_pL = (upperlimit_pL+lowerlimit_pL)/2.0;

% integral over theta_L, outer layer loop
neword_tL = conv1_tL*xg + conv2_tL;
mu_tL     = cos(neword_tL'); mu_tL = repmat(mu_tL,[ng 1]);
sin_tL    = sin(neword_tL'); sin_tL = repmat(sin_tL,[ng 1]);

% integral over phi_L, inner layer loop
neword_pL  = conv1_pL*xg + conv2_pL;
dotproduct1 = ( mu_tL*mu_tp + sin_tL*sin_tp.*cos(neword_pL-phiprime) );
dotproduct2 = ( mu_tL*mu_t  + sin_tL*sin_t .*cos(neword_pL-phi     ) );
dot_dotproduct12 = dotproduct1.*dotproduct2;
sum_pL = rho_Ld*wg.*hL/(2.0*pi).*abs(dot_dotproduct12);
sum_pL_tau = tau_Ld*wg.*hL/(2.0*pi).*abs(dot_dotproduct12);
sum_pL(dot_dotproduct12>0) = sum_pL_tau(dot_dotproduct12>0);
sum_pL = sum(sum_pL,1);
%end of inner layer loop

% finish the phi_L integral
sum_pL = sum_pL*conv1_pL;
sum_tL = sum(wg.*gL.*sum_pL');
% end of outer layer loop

% finish the theta_L integral
sum_tL   = sum_tL*conv1_tL;
Gamma_d  = sum_tL;

end

