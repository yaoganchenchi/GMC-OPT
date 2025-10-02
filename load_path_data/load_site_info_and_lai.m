function [site_info,site_name,LAI_this_site,LAI_cli,max_LAI_cli,max_month_cli,nyear_fluxdata]...
    =load_site_info_and_lai(fn,eeo,glc,isite)

site_info = loadMatData(fn);
site_name = site_info.ST_name{isite};
nyear_fluxdata = site_info.FULLSET_ED_YEAR(isite) - site_info.FULLSET_ST_YEAR(isite)+1;

% round LAI, with 0.05 intervals
site_info.LAI_stack=round(site_info.LAI_stack./eeo.DeltaL).*eeo.DeltaL;

% get climatology
LAI_this_site = site_info.LAI_stack(isite,:);
LAI_this_site = reshape(LAI_this_site,[glc.nmonth,length(LAI_this_site)/glc.nmonth]);
LAI_cli = mean(LAI_this_site,2,'omitnan');
[max_LAI_cli,max_month_cli]=max(LAI_cli);

end