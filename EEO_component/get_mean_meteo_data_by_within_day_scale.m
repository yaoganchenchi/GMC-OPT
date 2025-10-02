function mean_meteo_data=get_mean_meteo_data_by_within_day_scale...
    (within_day_scale,vis_hourly,sw_hourly,lw_hourly,PA_F,TA_F,SWC_F,...
    USTAR_F,ea_F,qa_F,glc,seb,eeo,day_mask)

% do calculation
valid_mask_hr = ~isnan(vis_hourly);
if strcmp(within_day_scale,'hourly')
    rhoa_F = (PA_F - ea_F)./(seb.R.*TA_F) + ea_F./(seb.Rv.*TA_F); % rho air

    Ta_hr = repmat(TA_F',[eeo.Nlayers,1]);Ta_hr(~valid_mask_hr) = nan;
    Pz_hr = repmat(PA_F',[eeo.Nlayers,1]);Pz_hr(~valid_mask_hr) = nan;
    qa_hr = repmat(qa_F',[eeo.Nlayers,1]);qa_hr(~valid_mask_hr) = nan;
    SWC_hr= repmat(SWC_F',[eeo.Nlayers,1]);SWC_hr(~valid_mask_hr) = nan;
    USTAR_hr = repmat(USTAR_F',[eeo.Nlayers,1]);USTAR_hr(~valid_mask_hr) = nan;
    rhoa_hr = repmat(rhoa_F',[eeo.Nlayers,1]);rhoa_hr(~valid_mask_hr) = nan;

    mean_meteo_data = struct('vis_ave',vis_hourly,'sw_ave',sw_hourly,'Pz_ave',Pz_hr,...
        'Ta_ave',Ta_hr,'SWC_ave',SWC_hr,...
        'qa_ave',qa_hr,'rhoa_ave',rhoa_hr,'USTAR_ave',USTAR_hr,...
        'lw_ave',lw_hourly);
end

end