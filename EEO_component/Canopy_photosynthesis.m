function [A,g_mol,gs,Vcmax,Jmax,J,Tleaf_opt,ci,lambda,WUE,Tr,PPFD,D,idx,ref_vcmax]=Canopy_photosynthesis(Ta,Tleaf,Pz,qa,SWC,ca,PAR,gb_c,opt_Vcmax,beta0,carbon,eeo)

%% start to calculate GPP
%% D
es = carbon.es_coef1.*exp(carbon.es_coef2.*(Tleaf-carbon.T_ref)./(Tleaf-carbon.T_ref+carbon.es_ref_temp2))./Pz;
ea = qa./carbon.gas_epsilon;
D = es-ea;
% Let VPD leaf to air > VPD air to air.
% 1. let leaf temp  = Ta
Ta2 = Ta.*ones(size(D));
es2 = carbon.es_coef1.*exp(carbon.es_coef2.*(Ta2-carbon.T_ref)./(Ta2-carbon.T_ref+carbon.es_ref_temp2))./Pz;
D2 = es2-ea;
D(D<D2) = D2(D<D2);
% 2. if it doesn't work, let D to be a tiny value greater than 0
D(D<0)=10^(-12);

%% Michaelis-Menten coefficient
Kc =carbon.Kc25.*Arrh(Tleaf,carbon.delHkc,carbon.T_25, carbon.R);  %exp(delHkc.*(T-T_ref) ./ (T_ref.*R.*T));
Ko =carbon.Ko25.*Arrh(Tleaf,carbon.delHko,carbon.T_25, carbon.R);  %exp(delHko.*(T-T_ref) ./ (T_ref.*R.*T));
K = Kc.*(1+carbon.oa./Ko);

%% ci, leaf CO2 concentration and lambda
Gamma_star = carbon.Gamma_star_25.*Arrh(Tleaf,carbon.delH_Gamma,carbon.T_25, carbon.R).*Pz./carbon.P0;
n_star = exp(-carbon.nB.*(Tleaf-carbon.T_25) ./(Tleaf+carbon.nC)./(carbon.T_25+carbon.nC));
inv_beta = carbon.inv_beta_prime.*n_star;

if eeo.linearity==1
    lambda = (1+exp(-beta0.*SWC))./(1-exp(-beta0.*SWC)).*ca.*inv_beta.*1./K ./ (1+sqrt(carbon.a.*inv_beta.*D./K)).^2;
    ci = ca.*(1-(carbon.a.*lambda.*D./ca).^0.5);
elseif eeo.linearity==0
    ci = ca./(1+(carbon.a.*inv_beta.*D.*carbon.n1./K).^0.5).*ones(size(SWC));
    lambda_temp=get_lambda_using_ci(carbon.a,D,K,ca,ci,Gamma_star);
    lambda = lambda_temp.l1.*(1+exp(-beta0.*SWC))./(1-exp(-beta0.*SWC)); % now, we use l1 solution, l2 is also acceptable and gives identical ci and gs
end
s=ci/ca; % added

PPFD = PAR.*carbon.c_PAR2PPFD;

phi_psii_500 = carbon.phi_psii_coef0+carbon.phi_psii_coef1 .* (Tleaf-carbon.T_ref) + carbon.phi_psii_coef2 .* (Tleaf-carbon.T_ref).^2; % 500 umol abs PPFD as a function of temperature
phi_psii_light = phi_psii_500.*(1-carbon.slope_photoinbit.*(PPFD./carbon.aL./eeo.cL-carbon.light_adapt_ref));
phi_psii_light(phi_psii_light<0)=10^-12;
phi=carbon.bL*phi_psii_light/carbon.hk; % must be less than 0.125 mol CO2 / mol Photon

%% Vcmax Jmax
% Vcmax type, 1: everymonth, 2: every year peak season,
if opt_Vcmax.acclimation_type==1 % monthly acclimation
    Vcmax_m = opt_Vcmax.Vcmax;
    Jmax_m = opt_Vcmax.Jmax;
    T_m = opt_Vcmax.Tleaf;
    
elseif opt_Vcmax.acclimation_type==2
    [nlayer_total,~] = size(PAR);
    nvalid_layer=sum(sum(PAR,2,'omitnan')>0);
    nvliad_max_layer = sum(~isnan(opt_Vcmax.Vcmax));
    Vcmax_m = nan(nlayer_total,1);
    Jmax_m = nan(nlayer_total,1);
    T_m = nan(nlayer_total,1);
    
    if opt_Vcmax.leaf_flush_type==1
        % assign Vcmax using uniform distribtuion according to the refrence Vcmax
        itv = round(linspace(1,nvliad_max_layer,nvalid_layer));
        temp=opt_Vcmax.Vcmax(itv);Vcmax_m(1:nvalid_layer) = temp(1:nvalid_layer);
        temp=opt_Vcmax.Jmax(itv);Jmax_m(1:nvalid_layer) = temp(1:nvalid_layer);
        temp=opt_Vcmax.Tleaf(itv);T_m(1:nvalid_layer) = temp(1:nvalid_layer);
        
    elseif opt_Vcmax.leaf_flush_type==2
        % assign Vcmax using from the bottom of the canopy
        Vcmax_m(1:nvalid_layer) = opt_Vcmax.Vcmax(nvliad_max_layer-nvalid_layer+1:nvliad_max_layer);
        Jmax_m(1:nvalid_layer) = opt_Vcmax.Jmax(nvliad_max_layer-nvalid_layer+1:nvliad_max_layer);
        T_m(1:nvalid_layer) = opt_Vcmax.Tleaf(nvliad_max_layer-nvalid_layer+1:nvliad_max_layer);
        
    elseif opt_Vcmax.leaf_flush_type==3
        % assign Vcmax using from the top of the canopy
        Vcmax_m = opt_Vcmax.Vcmax;
        Jmax_m = opt_Vcmax.Jmax;
        T_m = opt_Vcmax.Tleaf;
    end
end
ref_vcmax.Vcmax_m = Vcmax_m;
ref_vcmax.Jmax_m = Jmax_m;
ref_vcmax.T_m = T_m;

%% Vcmax
del_S = carbon.as-carbon.bs.*(Tleaf-carbon.T_ref); % intermedia variable
xx = Arrh(Tleaf,carbon.del_H,T_m, carbon.R).*(1+exp((T_m.*del_S-carbon.Hd)./(T_m.*carbon.R)))./(1+exp((Tleaf.*del_S-carbon.Hd)./(Tleaf.*carbon.R))); % intermedia variable
Vcmax = Vcmax_m.*xx;

%% Jmax and J
del_Sj = carbon.asj-carbon.bsj.*(Tleaf-carbon.T_ref);% intermedia variable
xxj = Arrh(Tleaf,carbon.del_Hj,T_m, carbon.R).*(1+exp((T_m.*del_Sj-carbon.Hd)./(T_m.*carbon.R)))./(1+exp((Tleaf.*del_Sj-carbon.Hd)./(Tleaf.*carbon.R)));% intermedia variable
Jmax = Jmax_m.*xxj;
J = carbon.hk.* phi.*PPFD ./(1+(4.* phi.*PPFD./Jmax).^2).^0.5;
%     J = (phi.*PPFD+Jmax - ((phi.*PPFD+Jmax).^2-4.*carbon.theta.*phi.*PPFD.*Jmax).^0.5)./(2.*carbon.theta);

%% stomatal conductance
g_c = carbon.n1.*Vcmax./(K+carbon.n1.*s.*ca).*(-1+(ca./(carbon.a.*lambda.*D)).^0.5);
g_j = carbon.n1.*J./carbon.hk./(carbon.n1.*s.*ca+2.*carbon.n2.*Gamma_star).*(-1+(ca./(carbon.a.*lambda.*D)).^0.5);

%% dark respiration
%Rd = coef_Rd.*Vcmax; %Following flux tower definition of GPP. Do not count Rd.

%% final photosynthesis
% Farquhar style
Ac = Vcmax.*(carbon.n1.*ci-carbon.n2.*Gamma_star) ./ (K+carbon.n1.*ci); % should not optimize for Rd
Aj = J./carbon.hk.*(carbon.n1.*ci-carbon.n2.*Gamma_star)./(carbon.n1.*ci+2.*carbon.n2.*Gamma_star);
Ac(Ac<0|D<0|ci<0)=0;
Aj(Aj<0|D<0|ci<0)=0;
temp=cat(3,Ac,Aj);
[A,idx] = min(temp,[],3,'omitnan');

valid_mask =  A>=0  & D>=0 ;%& Tleaf>T_ref;
A(~valid_mask)=nan; A(isnan(A) )=0; % set A=0 if pixel is veg but met condition is not suitable for photoysntheisis.
idx((isnan(A)|A==0))=nan;

g_mol=g_c;
g_mol(idx==2)=g_j(idx==2);

WUE = (ca.*lambda./carbon.a./D).^0.5; % umol C mol^-1 H2O
Tr = carbon.a.*g_mol.*D;
gs = g_mol.*(carbon.R.*Ta)./Pz; % Jones, 2014, p54 & Appendix 3
rs = 1./gs; % for CO2, (leaf boundary layer + stomatal) now in s/m, 
rs_photo = rs-1./gb_c; % remove the boundary layer resistance, 1.4 is CO2 and H2O ratio
gs = 1./rs_photo;% CO2
gs(gs<0) = 10^(-12); % must > 0
gs(~valid_mask)= 10^(-12); % if photosynthesis stops, transpiration should also stop
g_mol = gs ./ ((carbon.R.*Ta)./Pz); % stomatal conductance for CO2 after romoving leaf boundary layer resitance

Tleaf_opt = Tleaf; % Tleaf at n-1 th iteraction

end



