function GMCOPT_main(isite,iyear,imode,linearity,cali)
% This is the main functiion of GMCOPT
% Author: Chi Chen
% Date: 10/01/2025

if ~exist('isite','var')
    isite = 175; % when extended to global, this should be indexed by row and col number (or lat/long)
end
if ~exist('iyear','var')
    iyear = 2004;
end
if ~exist('imode','var')
    imode = 2; % mode 1, 2, 3, 4
end
if ~exist('linearity','var')
    linearity = 1;% 1: linear stomatal conductance model, 0: non-linear model
end
if ~exist('cali','var')
    cali = 3;  % 0: no cali; 1: bulk, 2: pft-specific, 3. site-specific
end

%%%% some model options
vis_type_dict = {'vis_median'};
within_day_scale_dict = {'hourly'}; % only hourly is avaliable in this version 
lad_dict={'Spherical','Planophile','Erectophile','Plagiophile','Extremophile','Uniform'};

vis_type = vis_type_dict{1}; % assume vis_median: 0.8 leaf absorption
within_day_scale = within_day_scale_dict{1}; % Vcmax acclimate to sub-daily scales
leaf_angle_distr = lad_dict{1}; % spherical leaf normal distribution

if imode==1
    acclimation_type = 1; % 1: monthly, 2: peak LAI month
    leaf_flush_type = 1;% 1: uniform, 2: keep botton, 3: keep top
elseif imode ==2
    acclimation_type = 2; 
    leaf_flush_type = 1;
    
elseif imode ==3
    acclimation_type = 2; 
    leaf_flush_type = 2;
elseif imode ==4
    acclimation_type = 2; 
    leaf_flush_type = 3;
end

fprintf('########################\n')
fprintf('INPUT PARAMETERS\n')
fprintf('iyear = %04d, isite = %03d\n',iyear,isite);
fprintf('cali = %d, linearity = %d, acclimation_type = %d (1: monthly, 2: peak LAI month), leaf_flush_type = %d (1: uniform, 2: keep botton, 3: keep top)\n',cali,linearity,acclimation_type,leaf_flush_type);
fprintf('vis_type = %s, within_day_scale = %s, leaf_angle_distr = %s\n',vis_type,within_day_scale, leaf_angle_distr);
fprintf('########################\n')

%% load paths and input data
A000_load_paths

% load some global constants
load_and_specify_constants
glc.iyear = iyear;
eeo.linearity = linearity;
fprintf('########################\n')
fprintf('Working Dir: %s\n',program_dir);
fprintf('########################\n')

%% Further loading paths and input data. In this section, ensure your input data is aligned with "this_year_data" and "glc" following the structure of the example inputs.
beta0=10^9; % currently not in use

% load site info, LAI, Lat, Lon
[site_info,site_name,LAI_this_site,LAI_cli,max_LAI_cli,max_month_cli,nyear_fluxdata]...
    =load_site_info_and_lai(site_info_path,eeo,glc,isite);
glc.LAI_this_site_this_year = LAI_this_site(:,glc.iyear-glc.MODIS_first_year+1);
glc.LON = site_info.LON(isite);
glc.LAT = site_info.LAT(isite);

% load calibration data
cali_data = loadMatData(sprintf('%scL_cali.mat',cali_root));
this_site_IGBP = site_info.IGBP{isite};
if cali==0 % Null model
    cL = 1.0;   
elseif cali==1 % bulk, max LAI month data
    cL = cali_data.hourly_maxLAI.bulk.slope;
elseif cali==2 % pft, max LAI month data
    cL = cali_data.hourly_maxLAI.(this_site_IGBP).slope;
    if isnan(cali_data.hourly_maxLAI.(this_site_IGBP).slope)
        cL = cali_data.hourly_maxLAI.bulk.slope;
    end
elseif cali==3 % site, max LAI month data
    cL = cali_data.hourly_maxLAI.site.slope(isite);
    if isnan(cali_data.hourly_maxLAI.site.slope(isite))
        cL = cali_data.hourly_maxLAI.(this_site_IGBP).slope;
    end
    if isnan(cali_data.hourly_maxLAI.(this_site_IGBP).slope)
        cL = cali_data.hourly_maxLAI.bulk.slope;       
    end  
end
if cL>1
    cL=1;
end
eeo.cL = cL;

% load EEO paths
load_EEO_paths

% load flux tower data
[flux_data_new,glc.nHour,glc.hh_intv,glc.nrow_per_year,glc.nrow_per_month,glc.syear,glc.eyear] =load_flux_data(meteo_root_dir,site_info,nyear_fluxdata,isite,site_name,glc);

% load global CO2 data
ca_NOAA = load_CO2(glc,CO2_path);
idx_this_year = (iyear-glc.MODIS_first_year)*glc.nmonth+1: (iyear-glc.MODIS_first_year+1)*glc.nmonth;
ca_this_year = ca_NOAA(idx_this_year);

% check whether the request year is within the available data years
if iyear < glc.syear || iyear > glc.eyear
    fprintf('Quit Program: The requested processing year (%04d) is outside the range of available data years (min:%04d,max:%04d). %d, %s,%s,%s\n',iyear,glc.syear,glc.eyear,isite,site_name,leaf_angle_distr,vis_type);
    return
end

% get the data of the requested year
idx_year = (iyear - glc.MODIS_first_year)*glc.nrow_per_year+1: (iyear - glc.MODIS_first_year+1)*glc.nrow_per_year;
fnames = fieldnames(flux_data_new);
this_year_data=[];
for ivar = 1:length(fnames)
    this_var = fnames{ivar};
    this_year_data.(this_var) = flux_data_new.(this_var)(idx_year,:);
end
if sum(~isnan(this_year_data.TA_F))==0 
    fprintf('Quit Program: No avaliable data to run the RT module (according to TA_F)\n');
    return
end

%% evaluate Red (PAR) and NIR radiation
fprintf('########################\n');
fprintf('Start processing the RT field.\n');
fprintf('For testing purposes in "A001_site_Rad_Field", it is recommended to set ng = 6 for the Gaussian quadrature.\n');
%%%% input vars: SW_IN, SW_IN_POT, PA, day_night_mask, LW_IN, Ts_Soil %%%%
%
RT_input_this_year_data.SW_IN_F = this_year_data.SW_IN_F;
RT_input_this_year_data.SW_IN_POT = this_year_data.SW_IN_POT;
RT_input_this_year_data.PA_F = this_year_data.PA_F;
RT_input_this_year_data.NIGHT = this_year_data.NIGHT;
RT_input_this_year_data.LW_IN_F=this_year_data.LW_IN_F;
RT_input_this_year_data.TS_F_MDS_1=this_year_data.TS_F_MDS_1; % soil surface temp
RT_input_this_year_data.TA_F = this_year_data.TA_F; % use it as leaf temperature (only for LW RT)

% RT for SW and LW
[SW_rad_field, LW_rad_field] = A001_site_Rad_Field(leaf_angle_distr,vis_type,eeo,glc,seb,RT_input_this_year_data);

rad_data = SW_rad_field;
rad_data.AB_lw_stack = LW_rad_field.AB_lw_stack;

fprintf('Finish processing the RT field.\n');
fprintf('########################\n');

%% run physiological and energy balance
fprintf('########################\n');
fprintf('Start processing the physiology and energy balance modules\n');

%%%% input variables (used): Ta, Pz, qa (or VPD), USTAR, ca, sw_abs, lw_abs, LAI %%%%
%%%% input variables (currently not necessary): SWC, LE, H, GPP_fluxnet (GPP_F) %%%%
this_year_data = rmfield(this_year_data,'TS_F_MDS_1');
this_year_data = rmfield(this_year_data,'LW_IN_F');
this_year_data = rmfield(this_year_data,'H_F_MDS');
this_year_data = rmfield(this_year_data,'SW_IN_POT');

% if cL = 1.0. All PAR_abs arrives chloroplast
data2save_hr=A002_site_Physio_EnergyBal_hourly(acclimation_type,leaf_flush_type,within_day_scale,carbon,seb,glc,eeo,rad_data,this_year_data,ca_this_year,beta0);
savename = sprintf('%s%s_hourly_%04d_mode%02d.mat',path_out,site_name,iyear,imode);
save(savename,'data2save_hr','-v7.3');

fprintf('Finish processing the physiology and energy balance modules\n');
fprintf('Check outputs: %s\n',path_out);
fprintf('########################\n');

end