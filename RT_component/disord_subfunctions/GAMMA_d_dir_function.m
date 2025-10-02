function Gamma_d_dir = GAMMA_d_dir_function (ng,xg,wg,gL,hL,rho_Ld,tau_Ld,muprime,phiprime,Gamma_d_dir)
% this routine evaluates the GAMMA_d function (muprime,phiprime ->
% mu,phi) where (mu,phi) are quadrature directions and checks for
% normalization

%------Variables ------
%   xg(ng), wg(ng)
% 	gL(ng), hL(ng)
%   muprime, phiprime, rho_Ld, tau_Ld
% 	upperlimit_p, lowerlimit_p, conv1_p, conv2_p
% 	mu, phi, Gamma_d_dir(ng,ng)

%  conversion factors to have the ordinates simulate phi
upperlimit_p = 2.0*pi;
lowerlimit_p = 0.0;
conv1_p = (upperlimit_p-lowerlimit_p)/2.0;
conv2_p = (upperlimit_p+lowerlimit_p)/2.0;

%  now get the Gamma_d_dir matrix direction by direction
for i = 1: ng
    mu = xg(i);
    for j = 1: ng
        phi = conv1_p*xg(j) + conv2_p;
        Gamma_d_dir(j,i)=GAMMA_d_function(ng,xg,wg,gL,hL,rho_Ld,tau_Ld,muprime,phiprime,mu,phi);
    end
end

end

