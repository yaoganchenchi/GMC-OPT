%% get yearly peak season data
function opt_Vcmax=Vcmax_yearly_peak_month(this_year_data,rad_data,...
    max_month_this_year,carbon,seb,glc,eeo,ca_this_year,beta0,within_day_scale)
%% load paths
A000_load_paths

%% default setting
opt_Vcmax.acclimation_type = 1; % monthly acclimation for this peak season calculation
opt_Vcmax.hr_daytime_mask = 1; %1: daytime (24-hour ave, daytime ave, noon)
imonth = max_month_this_year; % max LAI month in a year

%% get this month flux data
[glc,vis_hourly,sw_hourly,lw_hourly,...
    TA_F,SW_IN_F,PA_F,ea_F,qa_F,SWC_F,USTAR_F,ca,LE_F,...
    day_mask,GPP_F]...
    =get_this_month_fluxdata(this_year_data,rad_data,imonth,ca_this_year,glc,eeo);

nvalid_layer = sum(sum(vis_hourly,2,'omitnan')>0); % avaliable layers, determined by LAI of that month

%% get meteorology at required timescale (24-h, daytime, noon,hourly)
Vcmax_met=get_mean_meteo_data_by_within_day_scale...
    (within_day_scale,vis_hourly,sw_hourly,lw_hourly,...
    PA_F,TA_F,SWC_F,USTAR_F,ea_F,qa_F,glc,seb,eeo,day_mask);

%% currently muted (by setting beta0 with a very large value). In the future, this can be the leaf water potential inferred from soil moisture.
if strcmp(within_day_scale,'hourly')
    SWC_temp = SWC_F'.* linspace(glc.swc_min_pct,glc.swc_max_pct,nvalid_layer)';
    Vcmax_met.SWC_ave=nan(eeo.Nlayers,glc.nHour);Vcmax_met.SWC_ave(1:nvalid_layer,:)=SWC_temp;
else
    SWC_temp = Vcmax_met.SWC_ave.* linspace(glc.swc_min_pct,glc.swc_max_pct,nvalid_layer)';
    Vcmax_met.SWC_ave=nan(eeo.Nlayers,1);Vcmax_met.SWC_ave(1:nvalid_layer)=SWC_temp;
end

%% do an hourly calculation
%% convert rb heat to rb H2O vapor and CO2
Vcmax_met.rb_ave = 1/seb.Cv *(Vcmax_met.USTAR_ave/seb.dleaf).^(-0.5); % CLM5 technical note p60, for heat
Vcmax_met.rb_v_ave = Vcmax_met.rb_ave./(seb.Dv_Dh).^0.67; % Bonan 2014
Vcmax_met.gb_c_ave = 1./Vcmax_met.rb_v_ave ./seb.a_CO2_rb;

%% fill USTAR (in case USTAT at noon hour is missing)
if strcmp(within_day_scale,'noon') && isnan(Vcmax_met.USTAR_ave)
    USTAR_ave = mean(USTAR_F,'omitnan');
    Vcmax_met.rb_ave=1/seb.Cv *(USTAR_ave/seb.dleaf)^(-0.5);
    Vcmax_met.rb_v_ave = Vcmax_met.rb_ave./(seb.Dv_Dh)^0.67; % Bonan 2014
    Vcmax_met.gb_c_ave = 1./Vcmax_met.rb_v_ave ./seb.a_CO2_rb;
    fprintf('noon, USTAR filled for Vcmax opt calculation: year=%04d,imonth=%02d\n',glc.iyear,imonth)
end

% start energy balance iteration
if sum(Vcmax_met.SWC_ave(1,:),'omitnan')==0 || sum(Vcmax_met.Ta_ave(1,:),'omitnan')==0 ||...
        sum(Vcmax_met.Pz_ave(1,:),'omitnan')==0 || sum(Vcmax_met.qa_ave(1,:),'omitnan')==0 ...
        || sum(Vcmax_met.USTAR_ave(1,:),'omitnan')==0|| sum(Vcmax_met.lw_ave(1,:),'omitnan')==0
    fprintf('Vcmax calculation failed. Due to at least one of the followings are missing (SWC, Ta, Pz, qa, USTAR, LW_abs): year=%04d,imonth=%02d\n',glc.iyear,imonth)
    opt_Vcmax.valid=false;
else
    G=0; % assume no heat storage in leaves
    albedo=0; % already considered in sw input
    Tleaf = Vcmax_met.Ta_ave; % for the first guess, assume Tleaf=Ta
    
    for itr = 1: eeo.nitr
        % note that rs is for water vapor (gH2O = gCO2 * 1.6, rs = 1/gH2O)
        
        % physiology module
        [A,gs_CO2,Vcmax,Jmax,Tleaf_opt,PPFD_chl_opt,Tr]...
            = Canopy_photosynthesis_for_Vcmax(Vcmax_met.Ta_ave,Tleaf,...
            Vcmax_met.Pz_ave,Vcmax_met.qa_ave,Vcmax_met.SWC_ave,ca,...
            Vcmax_met.vis_ave/eeo.DeltaL*eeo.cL,Vcmax_met.gb_c_ave,...
            beta0,carbon,eeo);
        
        % energy balance module
        output=TRM_second(Vcmax_met.sw_ave./eeo.DeltaL, Vcmax_met.lw_ave./eeo.DeltaL./seb.emis_l, Vcmax_met.Ta_ave,...
            Vcmax_met.qa_ave, albedo, Vcmax_met.rb_ave,Vcmax_met.rb_v_ave,...
            1./(carbon.a.*gs_CO2), seb.emis_l, Vcmax_met.rhoa_ave, Vcmax_met.Pz_ave, G, seb);
        
        delta_Tleaf = output.Ts - Tleaf;
        
        %update Tleaf
        Tleaf = output.Ts;
        
        max_delta_Tleaf=max(abs(delta_Tleaf(:)),[],'omitnan');
        %fprintf('itr=%d,max delta_Tleaf = %f\n',itr,max_delta_Tleaf);
        if max_delta_Tleaf<eeo.eps_del_Tleaf && itr>1
            break
        end
    end
      
    No_A_flag = sum(A(:),'omitnan')==0 || isnan(sum(A(:),'omitnan'));   % determine if all A values are invalid or zero. If so, Vcmax.valid is false
    if No_A_flag
        opt_Vcmax.valid=false;
        
    elseif strcmp(within_day_scale,'hourly') && ~No_A_flag % subdaily hourly acclimation for vcmax
        Vcmax_out=nan(eeo.Nlayers,1);
        Jmax_out=nan(eeo.Nlayers,1);
        Tleaf_out = nan(eeo.Nlayers,1);
        PPFD_chl_out = nan(eeo.Nlayers,1);
        hour_idx = nan(eeo.Nlayers,1);
        
        for ilayer=1:nvalid_layer
            this_roi = A(ilayer,:);
            max_val = max(this_roi);
            if max_val==0
                continue % let all values to be nan
            end
            hour_idx(ilayer)=find(this_roi==max_val);
            Vcmax_out(ilayer)=Vcmax(ilayer,hour_idx(ilayer));
            Jmax_out(ilayer)=Jmax(ilayer,hour_idx(ilayer));
            Tleaf_out(ilayer)=Tleaf_opt(ilayer,hour_idx(ilayer));
            PPFD_chl_out(ilayer) = PPFD_chl_opt(ilayer,hour_idx(ilayer));
        end
        
        opt_Vcmax.Vcmax = Vcmax_out;
        opt_Vcmax.Jmax = Jmax_out;
        opt_Vcmax.Tleaf = Tleaf_out;
        opt_Vcmax.PPFD_chl = PPFD_chl_out;
        opt_Vcmax.hour_idx=hour_idx;
        opt_Vcmax.valid=true;
        opt_Vcmax.day_mask = day_mask;
    end
    
end
end