function Gamma_d_dif=GAMMA_d_dif_function(ng,xg,wg,gL,hL,rho_Ld,tau_Ld,Gamma_d_dif)
% this routine evaluates the Gamma_d_dif function for scattering
% from all quadrature directions to all exit quadrature directions

%------Variables ------
%   xg(ng), wg(ng), gL(ng), hL(ng), rho_Ld, tau_Ld
%   phiprime, muprime, Gdif(ng,ng)
%	upperlimit_pp, lowerlimit_pp, conv1_pp, conv2_pp
%	Gamma_d_dif(ng,ng,ng,ng), dummy(12,12)

Gamma_d_dir_dummy = nan(ng,ng);

% GAMMA_d_function conversion factors to have the ordinates simulate phiprime
upperlimit_pp = 2.0*pi;
lowerlimit_pp = 0.0;
conv1_pp = (upperlimit_pp-lowerlimit_pp)/2.0;
conv2_pp = (upperlimit_pp+lowerlimit_pp)/2.0;

%  now get the Gamma_d_dif matrix direction by direction
for i = 1: ng
    muprime = xg(i);
    for j = 1: ng
        phiprime = conv1_pp*xg(j) + conv2_pp;
        Gamma_d_dir_dummy = GAMMA_d_dir_function (ng,xg,wg,gL,hL,rho_Ld,tau_Ld,muprime,phiprime,Gamma_d_dir_dummy);
        
        for m = 1: ng
            for n = 1: ng
                Gamma_d_dif(j,i,n,m) = Gamma_d_dir_dummy(n,m);
            end
        end
    end
end
end

