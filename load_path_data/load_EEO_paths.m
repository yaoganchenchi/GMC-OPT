%% This script load paths that required by the EEO component

%% set paths and mkdir
if eeo.linearity==1
    str = 'linear';
elseif eeo.linearity==0
    str = 'nonlinear';
end

if cali==0
    cali_str = '100';
elseif cali==1
    cali_str = 'bulk_maxLAI';
elseif cali==2
    cali_str = 'pft_maxLAI';
elseif cali==3
    cali_str = 'site_maxLAI';
end

path_out = sprintf('%sOutput/%s_cL%s/%s_%s/diel_monthly_GPP_hourly/%s/',program_dir,str,cali_str,vis_type,within_day_scale,leaf_angle_distr);
if ~(exist(path_out,'dir') ==7)
    mkdir(path_out);
end