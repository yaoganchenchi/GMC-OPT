%% declare output variables
function [data2save_hr] = declare_model_output_variable(Nlayers,nHour,nmonth)

%% declare outputs for leaf level
% for global model output
data2save_hr.leaf.A = nan(Nlayers,nHour,nmonth); % μmol CO2 m^-3 s^-1(per unit LAI)
data2save_hr.leaf.Tr = nan(Nlayers,nHour,nmonth); % mol H2O m^-3 s^-1 (per unit LAI)
data2save_hr.leaf.Vcmax_25 = nan(Nlayers,nmonth); % at 25 deg, μmol C m^-3 s^-1 (per unit LAI)
data2save_hr.leaf.Jmax_25 = nan(Nlayers,nmonth); % at 25 deg, μmol electron m^-3 s^-1 (per unit LAI)
data2save_hr.leaf.gsCO2_mol = nan(Nlayers,nHour,nmonth);% two-side leaf, mol CO2 m^-3 s^-1 (per unit LAI)
data2save_hr.leaf.ci = nan(Nlayers,nHour,nmonth); % μmol CO2 mol^-1 air
data2save_hr.leaf.lambda = nan(Nlayers,nHour,nmonth); % μmol CO2 mol^-1 H2O, (~on the order of 10^3)
data2save_hr.leaf.PPFD_abs_chl = nan(Nlayers,nHour,nmonth); % mol photon m^-2 s^-1, to get W m^-3 (per unit LAI),divide by 4.55

data2save_hr.leaf.Tleaf4A=nan(Nlayers,nHour,nmonth); % K
data2save_hr.leaf.H = nan(Nlayers,nHour,nmonth); % W m^-3 (per unit LAI)
data2save_hr.leaf.LE = nan(Nlayers,nHour,nmonth);% W m^-3 (per unit LAI)

data2save_hr.leaf.PAR_abs = nan(Nlayers,nHour,nmonth); % W m^-3 (per unit LAI)
data2save_hr.leaf.SW_abs = nan(Nlayers,nHour,nmonth); % W m^-3 (per unit LAI)
data2save_hr.leaf.LW_abs = nan(Nlayers,nHour,nmonth); % W m^-3 (per unit LAI)

data2save_hr.leaf.Vcmax=nan(Nlayers,nHour,nmonth); % at operating temperature
data2save_hr.leaf.Jmax=nan(Nlayers,nHour,nmonth); % at operating temperature
data2save_hr.leaf.J=nan(Nlayers,nHour,nmonth); % at operating temperature
data2save_hr.leaf.Vcmax_m = nan(Nlayers,nmonth); % at reference optimal temperature
data2save_hr.leaf.Jmax_m = nan(Nlayers,nmonth); % at reference optimal temperature
data2save_hr.leaf.T_m = nan(Nlayers,nmonth); % the reference optimal temperature
data2save_hr.leaf.hour_idx = nan(Nlayers,nmonth); % hour idx

data2save_hr.leaf.Rn = nan(Nlayers,nHour,nmonth); %(per unit LAI)
data2save_hr.leaf.energy_balance=nan(Nlayers,nHour,nmonth); %(per unit LAI)

%% declare outputs for canopy level
% for global model output
data2save_hr.cano.GPP = nan(nHour,nmonth); % μmol CO2 m^-2 s^-1
data2save_hr.cano.Tr = nan(nHour,nmonth); % mol H2O m^-2 s^-1
data2save_hr.cano.Vcmax_25_wA = nan(nHour,nmonth);% μmol CO2 m^-2 s^-1
data2save_hr.cano.Jmax_25_wA = nan(nHour,nmonth);%μmol electron m^-2 s^-1
data2save_hr.cano.gsCO2_mol = nan(nHour,nmonth);% two-side leaf, mol CO2 m^-2 s^-1

data2save_hr.cano.FPAR = nan(nHour,nmonth); % unitless

% for internal analysis use
data2save_hr.cano.Vcmax=nan(nHour,nmonth); % μmol CO2 m^-2 s^-1
data2save_hr.cano.Jmax=nan(nHour,nmonth); %μmol electron m^-2 s^-1
data2save_hr.cano.J=nan(nHour,nmonth);%μmol electron m^-2 s^-1
data2save_hr.cano.Vcmax_wA=nan(nHour,nmonth); %weighted by GPP, % μmol CO2 m^-2 s^-1
data2save_hr.cano.Jmax_wA=nan(nHour,nmonth);%weighted by GPP, μmol electron m^-2 s^-1
data2save_hr.cano.J_wA=nan(nHour,nmonth);%weighted by GPP, μmol electron m^-2 s^-1

data2save_hr.cano.Vcmax_25 = nan(nHour,nmonth); % μmol CO2 m^-2 s^-1
data2save_hr.cano.Jmax_25 = nan(nHour,nmonth); %μmol electron m^-2 s^-1

%% declare outputs for input data at canopy level
data2save_hr.inputs.DeltaL = nan;
data2save_hr.inputs.fchl = nan;
data2save_hr.inputs.LAI=nan(1,nmonth);
data2save_hr.inputs.maxLAImonth=nan(1,nmonth);
data2save_hr.inputs.maxLAI=nan(1,nmonth);

data2save_hr.inputs.SW=nan(nHour,nmonth);
data2save_hr.inputs.Ta=nan(nHour,nmonth);
data2save_hr.inputs.Pz=nan(nHour,nmonth);
data2save_hr.inputs.qa=nan(nHour,nmonth);
data2save_hr.inputs.SWC=nan(Nlayers,nHour,nmonth);
data2save_hr.inputs.ca=nan(1,nmonth);

data2save_hr.inputs.GPP_flux=nan(nHour,nmonth);

%%%%%%
end