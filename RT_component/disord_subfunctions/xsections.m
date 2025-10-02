%% xsection sub functions
function [xg,wg,Gdir,Gdif,Gamma_d_dir,Gamma_d_dif]=xsections(ng,muprime,phiprime,rho_Ld,tau_Ld,lnd)
% this simple program illustrates

% (1) obtain quadrature [ng,xg,wg]
% (2) obtain leaf normal orientation pdfs [gL,hL]
% (3) evaluate the G function [G(mu,phi)]
% (4) evaluate the Gamma_d function [Gamma_d(phip,mup->phi,mu]

% begin declarations
% Gdir: a single value
%xg = nan(ng,1); wg = nan(ng,1);
gL = nan(ng,1); hL = nan(ng,1);
Gdif = nan(ng,ng);
Gamma_d_dir = nan(ng,ng);
Gamma_d_dif = nan(ng,ng,ng,ng);
% end declarations

%% get quadrature
%[xg,wg]=gauss_quad(ng,xg,wg);
[xg,wg]=gauss_quad_new(ng);

%% get pdf of leaf normal orientation gL(thetaL) and hL(phiL)
[gL,hL] = leaf_normal_pdf(xg,wg,gL,hL,lnd);

%% get G_FUNCTION for a direction OMEGA^prime
Gdir=G_dir_function(ng,xg,wg,gL,hL,muprime,phiprime);

%% get G_FUNCTION for all quadrature directions
% Gdif(phi,mu)
Gdif=G_dif_function(ng,xg,wg,gL,hL,Gdif);

%% get GAMMA_d (phiprime, muprime -> phi, mu) where (mu,phi) are all quadrature directions
%  Gamma_d_dir(phi,mu)
Gamma_d_dir = GAMMA_d_dir_function (ng,xg,wg,gL,hL,rho_Ld,tau_Ld,muprime,phiprime,Gamma_d_dir);

%% get GAMMA_d (muprime,phiprime -> mu,phi) where both the incident and exit directions are all quadrature directions
%Gamma_d_dif(phiprime,muprime,phi,mu)
Gamma_d_dif=GAMMA_d_dif_function(ng,xg,wg,gL,hL,rho_Ld,tau_Ld,Gamma_d_dif);

end

