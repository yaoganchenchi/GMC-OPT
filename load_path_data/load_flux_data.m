function [flux_data_new,nHour,hh_intv,nrow_per_year,nrow_per_month,syear,eyear]=load_flux_data(meteo_root_dir,site_info,nyear_fluxdata,isite,site_name,glc)
% syear:  the first year of your data or MODIS data, wihch ever is later
% eyear: the last year of your data or MODIS data, which ever is earlier

fn = sprintf('%s%s_MM_DayTime_%04d_%04d.mat',meteo_root_dir,site_name,site_info.FULLSET_ST_YEAR(isite),site_info.FULLSET_ED_YEAR(isite));
flux_data = loadMatData(fn);

% GPP_F is used as a reference. It represents the average of all available FluxNet GPP estimates. 
% You may provide one GPP input solely for reference purposes.
flux_data.GPP_F = mean([flux_data.GPP_DT_CUT_MEAN,flux_data.GPP_DT_VUT_MEAN,flux_data.GPP_NT_CUT_MEAN,flux_data.GPP_NT_VUT_MEAN],2);

%fnames = fieldnames(flux_data);
fnames = {'TA_F', 'SW_IN_F', 'SW_IN_POT', 'PA_F', 'VPD_F', 'SWC_F_MDS_1', ...
    'LE_F_MDS', 'H_F_MDS', 'TS_F_MDS_1', 'LW_IN_F', 'USTAR', 'NIGHT','GPP_F'};
nrow = length(flux_data.(fnames{1}));

nHour = nrow / glc.nmonth / nyear_fluxdata;
hh_intv = 24 / nHour; % one day,24 hours
nrow_per_year = nHour*glc.nmonth;
nrow_per_month = nHour;

if glc.MODIS_first_year>=site_info.FULLSET_ST_YEAR(isite)
    syear = glc.MODIS_first_year;
    eyear = site_info.FULLSET_ED_YEAR(isite);
    idx_use1 = (glc.MODIS_first_year - site_info.FULLSET_ST_YEAR(isite))*nrow_per_year +1: (site_info.FULLSET_ED_YEAR(isite) - site_info.FULLSET_ST_YEAR(isite)+1)*nrow_per_year;
    idx_use2= 1:length(idx_use1);
elseif glc.MODIS_first_year<site_info.FULLSET_ST_YEAR(isite)
    syear = site_info.FULLSET_ST_YEAR(isite);
    eyear = site_info.FULLSET_ED_YEAR(isite);
    idx_use1 = 1: (site_info.FULLSET_ED_YEAR(isite) - site_info.FULLSET_ST_YEAR(isite)+1)*nrow_per_year;
    idx_use2= (site_info.FULLSET_ST_YEAR(isite)-glc.MODIS_first_year)*nrow_per_year +1:(site_info.FULLSET_ST_YEAR(isite)-glc.MODIS_first_year)*nrow_per_year+length(idx_use1);
end

flux_data_new=[];
for ivar = 1:length(fnames)
    this_var = fnames{ivar};
    flux_data_new.(this_var) = nan(glc.nyear*nrow_per_year,1);
    flux_data_new.(this_var)(idx_use2,:) = flux_data.(this_var)(idx_use1,:);
end

end