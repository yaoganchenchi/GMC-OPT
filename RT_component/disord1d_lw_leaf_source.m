function AB_nlayer=disord1d_lw_leaf_source(theta_o,fdir,LAI,lnd,wavelength,ng)
% longwave radiation field part 2. Jump start to compute leaf as a emission
% source

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
if ~exist('wavelength','var')
    wavelength = 'lw_down';
end

%   standard constants
degtorad = pi/180.0;
epsilon = 0.0001;

%BEGIN INPUTS
phi_o   = 0.0;
theta_o = theta_o*degtorad;
phi_o   = phi_o*degtorad;
mu_o    = cos(theta_o);
Ftot    = 1.0;
DeltaL  = 0.05;
I_d     = Ftot*(1-fdir)/pi/4;% leaf as partical source
Nlayers = round(LAI/DeltaL);

if strcmp(wavelength,'lw_down')
    rho_Ld  = 0.01;
    tau_Ld  = 0.01;
    R_s     = 0.04;   
elseif strcmp(wavelength,'lw_up')
    rho_Ld  = 0.01;
    tau_Ld  = 0.01;
    R_s     = 0.00;
elseif strcmp(wavelength,'leaf_source')
    rho_Ld  = 0.01;
    tau_Ld  = 0.01;
    R_s     = 0.04;
end
% END INPUTS

I_o_uc_d = zeros(Nlayers,1); % keep zero
I_o_uc_u=zeros(Nlayers,ng,ng); % keep zero

I_d_uc_d=zeros(Nlayers,ng,ng); % keep zero
I_d_uc_u=zeros(Nlayers,ng,ng);  % keep zero    
S=zeros(Nlayers,ng,ng);
Ic=zeros(Nlayers+1,ng,ng);

%% get cross sections
[xg,wg,Gdir,Gdif,Gamma_d_dir,Gamma_d_dif]=xsections(ng,mu_o,phi_o,rho_Ld,tau_Ld,lnd);

%% Evaluate First-Collision Source Q
Q = I_d * ones(Nlayers,ng,ng); % 4pi whole sphere

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