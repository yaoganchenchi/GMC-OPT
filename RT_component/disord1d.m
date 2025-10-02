function AB_nlayer=disord1d(theta_o,fdir,LAI,lnd,vis_nir,ng)
% Numerical solution of the 1D RTE with the discrete ordinate finite kernel method
% Myneni, Ranga B., Ross, J., & Asrar, G. (1989). A review on the theory of photon transport in leaf canopies. Agricultural and Forest Meteorology, 45(1–2), 1–153. 

if ~exist('theta_o','var')
    theta_o = 120;
end
if ~exist('fdir','var')
    fdir = 0.7;
end
if ~exist('LAI','var')
    LAI = 3.0;
end
if ~exist('lnd','var')
    lnd = 'Planophile';
end
if ~exist('vis_nir','var')
    vis_nir = 'vis_median';
end

%   standard constants
degtorad = pi/180.0;
epsilon = 0.0001;

% ___________________________________________________________________
% ng            : gauss quadrature order (4,6,8,10 or 12)
% Nlayers       : number of spatial nodes
% Epsilon       : convergence criterion for iteration of the
%                 scattering integral
% xg            : gauss ordinates (between -1 to 1)
% wg            : gauss weights (between 0 and 1)
% gL            : pdf of leaf normal inclination (assumed planophile)
% hL            : pdf of lear normal azimuths (assumed uniform)
% theta_o       : polar angle of the sun (between 90 and 180 degrees)
% phi_o         : azimuthal angle of the sun (between 0 and 360
%                  degrees)
% mu_o          : cosine of theta_o
% Ftot          : total incident flux density (1 W/m2)
% fdir          : fraction of Ftot in the direct solar beam
% I_o           : intensity of the direct beam
% I_d           : intensity of diffuse sky light (assumed isotropic)
% rho_Ld        : leaf hemispherical reflectance (diffuse internal
%                 scattering)
% tau_Ld        : leaf hemispherical transmittance (diffuse internal
%                 scattering)
% LAI           : one-sided leaf area per unit ground area
% DeltaL        : thickness of sptatial cells (LAI/Nlayers)
% R_s           : soil hemipsherical reflectance (assumed Lambertian)
% Gdir          : Geometry factor for direct solar radiation
% Gdif          : Geometry factor for scattered radiation field
% dummy         : well, a dummy variable
% Gamma_d_dir   : area scattering phase function for direct solar
%                 radiation
% Gamma_d_dif   : area scattering phase function for scattered
%                 radiation field
% I_o_uc_d      : downward uncollided direct solar radiation
% I_o_uc_u      : upward uncollided direct solar radiation
% I_d_uc_d      : downward uncollided diffuse sky radiation
% I_d_uc_u      : upward uncollided diffuse sky radiation
% F_o_uc_d_soil : donward uncollided flux density of direct solar
%                 radiation incident on the ground below the canopy
% F_d_uc_d_soil : donward uncollided flux density of diffuse sky
%                 radiation incident on the ground below the canopy
% Q             : first collision source
% S             : multiple collision source
% Ic            : collided intensity field
% Ic_old        : canopy leaving collided intensity field from
%                 previous iteration
% AB_nlayer     : absorptance of each layer
% convergence   : logical flag to test for convergence of the
%                 iteration on the multiple collision source
% theta_v       : view polar angle
% phi_v         : view azimuthal angle
% ___________________________________________________________________

%BEGIN INPUTS
%theta_o = 120.0;
phi_o   = 0.0;
theta_o = theta_o*degtorad;
phi_o   = phi_o*degtorad;
mu_o    = cos(theta_o);
Ftot    = 1.0;
%fdir    = 0.7;
I_o     = Ftot*(fdir/(abs(mu_o)));
I_d     = Ftot*(1-fdir)/pi ;
%LAI     = 3.0;
DeltaL  = 0.05;
Nlayers = round(LAI/DeltaL);
%lnd = 'Planophile';

if strcmp(vis_nir,'nir')
    % NIR
    rho_Ld  = 0.475;
    tau_Ld  = 0.45;
    R_s     = 0.2;
elseif strcmp(vis_nir,'vis_median')
    % Red
    rho_Ld  = 0.1;
    tau_Ld  = 0.1;
    R_s     = 0.125;
end
% END INPUTS

%I_o_uc_d = zeros(Nlayers,1);
%Q=zeros(Nlayers,ng,ng);  
I_o_uc_u=zeros(Nlayers,ng,ng);
I_d_uc_d=zeros(Nlayers,ng,ng); 
I_d_uc_u=zeros(Nlayers,ng,ng);      
S=zeros(Nlayers,ng,ng);
Ic=zeros(Nlayers+1,ng,ng);

%% get cross sections
[xg,wg,Gdir,Gdif,Gamma_d_dir,Gamma_d_dif]=xsections(ng,mu_o,phi_o,rho_Ld,tau_Ld,lnd);

%% Evaluate Uncollided direct solar radiation (I_o_uc_d, I_o_uc_u)
I_o_uc_d=I_o_uncol_down(DeltaL,Nlayers,Gdir,I_o,mu_o);
I_o_uc_u=I_o_uncol_up(DeltaL,Nlayers,Gdir, Gdif,I_o,mu_o,ng,xg,R_s,I_o_uc_u);% check if F_o_uc_d_soil is needed

%% Evaluate Uncollided diffuse sky radiation (I_d_uc_d, I_d_uc_u)
I_d_uc_d=I_d_uncol_down(DeltaL,Nlayers,Gdif,I_d,ng,xg,I_d_uc_d);
I_d_uc_u=I_d_uncol_up(DeltaL,Nlayers,Gdif,I_d,ng,xg,wg,R_s,I_d_uc_u); % check if F_d_uc_d_soil is needed

%% Evaluate First-Collision Source Q
Q=FCS(Nlayers,ng,wg,Gamma_d_dir,Gamma_d_dif,I_o_uc_d, I_o_uc_u,I_d_uc_d, I_d_uc_u);

%% Iterate on the Multiple-Collision Source S
for ims = 1: 100
    % add multiple-collision source to first-collision source
    S = Q+S;
    
    % sweep downwards plus handle the bottom boundary condition
    Ic=SWEEP_DOWN(Nlayers,ng,xg,wg,Gdif,DeltaL,S,R_s,Ic);
    
    % sweep upwards and check for convergence
    [Ic,convergence]=SWEEP_UP(Nlayers,ng,xg,Gdif,DeltaL,S,Ic,epsilon);
    
    % evaluate Multiple-Collision Source
    if (convergence)
        break;
    end
    [S]=MULTI_COLL_S(Nlayers,ng,wg,Gamma_d_dif,Ic,S);
    %fprintf('ims = %d\n',ims)
end

%% do each layer
[~,~,AB_nlayer]=ENERGY_BAL_Nlayer(Nlayers,ng,wg,DeltaL,rho_Ld,tau_Ld,Gdir,Gdif,I_o_uc_d,I_d_uc_d,I_o_uc_u,I_d_uc_u,Ic);

end

% Keep this line