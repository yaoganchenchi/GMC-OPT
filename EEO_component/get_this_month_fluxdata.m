function [glc,vis_hourly,sw_hourly,lw_hourly,...
    TA_F,SW_IN_F,PA_F,ea_F,qa_F,SWC_F,USTAR_F,ca,LE_F,...
    day_mask,GPP_F]=get_this_month_fluxdata(this_year_data,rad_data,imonth,ca_this_year,glc,eeo)

idx_month = (imonth-1)*glc.nrow_per_month+1:imonth*glc.nrow_per_month;
this_month_data=[];
fnames = fieldnames(this_year_data);
for ivar = 1:length(fnames)
    this_var = fnames{ivar};
    this_month_data.(this_var) = this_year_data.(this_var)(idx_month,:);
end


TA_F = this_month_data.TA_F+273.15;
SW_IN_F = this_month_data.SW_IN_F;
PA_F = this_month_data.PA_F*1000; %kPa
VPD_F = this_month_data.VPD_F; % hPa
%SWC_F = this_month_data.SWC_F_MDS_1/100;
SWC_F = ones(size(TA_F)); % currently not in use. need to switch once plant hydraulics is used


LE_F = this_month_data.LE_F_MDS; % Wm^-2
USTAR_F = this_month_data.USTAR;
ca = ca_this_year(imonth); % share the same idx as LAI, 2001-2014
GPP_F = this_month_data.GPP_F;% merged in load_flux_data

esat_F = es(TA_F);
ea_F = esat_F - VPD_F.*100;
qa_F = ea_F./PA_F.*eeo.gas_epsilon;
qa_F(qa_F<0)=nan;

% get daytime mask, at hourly or half-hourly frequency
vis_hourly = rad_data.AB_vis_I_stack(:,:,imonth);
nir_hourly = rad_data.AB_nir_I_stack(:,:,imonth);
sw_hourly = vis_hourly +nir_hourly;
lw_hourly = rad_data.AB_lw_stack(:,:,imonth);
day_mask = ~isnan(vis_hourly(1,:))';%day time mask
glc.day_length_ratio = sum(day_mask)/length(day_mask);% day_length ratio

end