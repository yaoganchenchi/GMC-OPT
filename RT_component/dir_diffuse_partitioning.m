function [fdir_vis,fdir_nir]=dir_diffuse_partitioning(SW_IN_F,SW_IN_POT,PA_F,SZA_ave)
% A., W. & J.M., N. Partitioning solar radiation into direct and diffuse, 
% visible and near-infrared components. Agr Forest Meteorol 34, 205–213 (1985).
% SZA must be less than 80 deg

A=0.9;
B=0.7;
C=0.88;
D=0.68;
P0 = 101.325;
RDV = 600.*cos(SZA_ave).*exp(-0.185.*PA_F/P0./cos(SZA_ave)); % direct visible beam radiation
RdV = 0.4*(600-RDV).*cos(SZA_ave); % Potential visible diffuse radiation

xxx = -1.1950+0.4459.*log10(1./cos(SZA_ave))-0.0345.*log10(1./cos(SZA_ave)).^2;
omega = 1320.* 10.^xxx;
RDN = (720.*exp(-0.06.*PA_F/P0./cos(SZA_ave))-omega).*cos(SZA_ave); % direct NIR beam radiation
RdN = 0.6*(720-RDN-omega).*cos(SZA_ave);% Potential NIR diffuse radiation

RATIO = SW_IN_F./SW_IN_POT;

% fraction of direct visible radiation
fdir_vis = RDV./(RDV+RdV).*(1-((A-RATIO)./B).^2/3);

% fraction of direct NIR radation
fdir_nir = RDN./(RDN+RdN).*(1-((C-RATIO)./D).^2/3);

% exclude some extreme case, when SZA is high.
fdir_vis(fdir_vis<0)=0;
fdir_nir(fdir_nir<0)=0;

fdir_vis(fdir_vis>1)=1;
fdir_nir(fdir_nir>1)=1;
end