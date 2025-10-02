function [carbon,seb,glc,eeo]= getConstants()

%% energy balance constants
seb.cp  = 1004.64;      % specific heat at constant pressure, J/kg/K
seb.Lv  = 2.501*10^6; % latent heat of vaparization
seb.R   = 287.058;       % dry air gas constant J/kg/K
seb.Rv  = 461.5;       % water vapor gas constant J/kg/K
seb.g   = 9.80616;      % gravity constant m/s^2
seb.k   = 0.41 ;      % von-Karman constant
seb.emis_l  = 0.98 ; % emissivity, Determining the emissivity of the leaves of nine horticultural crops by means of infrared thermography
seb.emis_s  = 0.96 ; % soil emissivity, CLM Tech Note p41
seb.sb          = 5.670367*10^(-8) ; % stephan-boltzman constant, W/(m^2 K^4)
seb.epsilon= seb.R/seb.Rv  ; 

% constants for leaf boundary layer resistance
seb.Cv = 0.01; % two side, m s^-0.5, CLM5 Technical Note p60, or should be cv =1/200 if one-side leaf 
seb.dleaf = 0.04; % m, CLM5 Technical Note p61
seb.Dv_Dh = 1.15; % from Bonan 2014
seb.a_CO2_rb = 1.4;%        relative diffusivity of water vapour w.r.t carbon dioxide for boundary layer resitance      

%% constant for carbon
carbon.gas_epsilon = seb.epsilon;       % gas_epsilon    a constant                                      0.622
carbon.es_coef1 = 611;            % es_coef1       coefficient for saturation vapor pressure       611              Eq. 3.9a from Dingman in Pa
carbon.es_coef2 = 17.27;          % es_coef2       coefficient for saturation vapor pressure       17.27
carbon.T_ref = 273.15;            % T_ref          temperature for freezing point of water         273.15 K
carbon.es_ref_temp2 = 237.3;      % es_ref_temp2   coefficient for saturation vapor pressure       237.3 K

carbon.T_25 = 298.15; % for Vcmax reference temperature
carbon.del_H = 71513;   % activation energy for vcmax
carbon.as = 668.39; % intercept for vcmax
carbon.bs=1.07; % slope for vcmax
carbon.Hd = 200000;  % deactivation energy for vcmax/jmax
carbon.del_Hj = 49884; % for Jmax
carbon.asj= 659.7; % intercept for Jmax
carbon.bsj = 0.75; % slope for Jmax
carbon.c_star = 0.41;% Wang et al 2017

carbon.f_vis = 0.45;% f_vis          fraction of visible radiation in SW             0.45
carbon.c_PAR2PPFD = 4.55;% c_PAR2PPFD     coefficient to convert PAR to PPFD              4.55 umol J^-1
carbon.P0 = 101325; % P0             pressure at sea level                           101325 Pa, if needed
carbon.Kc25 = 404.9; % Kc25           Kc at 25 C                                      404.9*10^-6*P0          Bernacchi et al 2001, Stocker et al. (2020)
carbon.Ko25 = 278.4; % Ko25           Ko at 25 C                                      278.4*10^-3*P0          
carbon.delHkc = 79430; % delHkc         activation energy for Kc in Arrhenius equation  79430 J mol^-1    
carbon.delHko = 36380; % delHko         activation energy for Ko in Arrhenius equation  36380 J mol^-1    

carbon.R = 8.3145; % R              universal gas constant                          8.3145 J mol^-1 K^-1
carbon.inv_beta_prime = 1/146*10^6; % inv_beta       Prentice's b/a                                  1/146   Stocker et al. (2020)
carbon.nB=580; %vogel equation
carbon.nC=-138;

carbon.a= 1.6;% a              relative diffusivity of water vapour w.r.t carbon dioxide for stomatal conductance    1.6
carbon.s = 0.7; % s            initial long-term average ci/ca ratio (you can randomly assign one, just for place holding, this will be re-calculated by the model anyway)
carbon.oa = 209.5; % oa             oxygen concentration                         mmol mol^-1
carbon.delH_Gamma = 37830; % delH_Gamma     activation energy for Kc in Arrhenius equation  37830 J mol^-1     Bernacchi et al 2001
carbon.Gamma_star_25 = 42.75;% Gamma_star_25  CO2 compensation point at 25C                   42.75 umol mol^-1  Bernacchi et al 2001

carbon.aL = 0.80; % = 1-rho-tau, according to leaf single scattering albedo (for a typical leaf, the medium scenerio in this RT experiment)
carbon.hk = 4; % four electrons to reduce one CO2 molecule 
carbon.bL = 0.5; % fraction of absorbed photos reached chloraplast that reaches PSII
carbon.phi_psii_coef0 = 0.352;
carbon.phi_psii_coef1 = 0.022;
carbon.phi_psii_coef2 = -0.00034; % Bernacchi et al 2003
carbon.phi_psii_max_dark = 0.85;
carbon.slope_photoinbit = 1/12000 * (carbon.hk / carbon.bL / carbon.aL); % Long 1994; Werner 1996; 1/12000 * (glc.hk / glc.bL / glc.aL)
carbon.light_adapt_ref = 500; %umol w m^-2 s^-1

carbon.coef_Rd = 0.015;%*ones(size(LC)); % coef_Rd        0.015 for C3, 0.015 for C4                      Collatz et al 1991
carbon.umol2gc = 12.0107*10^-6; % umol CO2 to gram Carbon
carbon.mol2gH20=18.0153; % mol H2O to gram H2O

%% for the eeo model iteration
eeo.nitr = 100; % max iteration times
eeo.eps_del_Tleaf = 0.1; % minimum leaf T error, K
eeo.DeltaL  = 0.05; % canopy layer interval
eeo.Nlayers = 7/eeo.DeltaL; %LAI 7 is max
eeo.gas_epsilon=seb.epsilon;

%% some global constants
glc.day_month = [31,28,31,30,31,30,31,31,30,31,30,31]; % no leap year considered
glc.sec_day=86400;
glc.nsecond_per_hour = 3600;
glc.nmonth=12;
end


