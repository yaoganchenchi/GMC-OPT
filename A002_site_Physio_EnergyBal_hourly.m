function data2save_hr=A002_site_Physio_EnergyBal_hourly(acclimation_type,leaf_flush_type,within_day_scale,carbon,seb,glc,eeo,rad_data,this_year_data,ca_this_year,beta0)
% Aug-19-2025

if ~exist('acclimation_type','var')
    acclimation_type = 2; %1: monthly, 2: yearly peak month
end
if ~exist('leaf_flush_type','var')
    leaf_flush_type = 1;
end

if ~exist('within_day_scale','var')
    within_day_scale = 'hourly';
end

% load paths and add dirctory for model tools
A000_load_paths

% declare saving outputs for daytime average and hourly
[data2save_hr] = declare_model_output_variable(eeo.Nlayers,glc.nHour,glc.nmonth);

%% do monthly calculation
% get the max LAI month of this year
[max_LAI_this_year,max_month_this_year]=max(glc.LAI_this_site_this_year);

% Opt Vcmax based on the max LAI month, only needed when using peak month acclimation approach
opt_Vcmax=Vcmax_yearly_peak_month(this_year_data,rad_data,...
    max_month_this_year,carbon,seb,glc,eeo,ca_this_year,beta0,within_day_scale);
opt_Vcmax.leaf_flush_type=leaf_flush_type;

if  ~opt_Vcmax.valid
    % valid=0, skip; valid =1, keep and not skip
    fprintf('Quit Program: opt_Vcmax is not avaliable for this year (%04d)\n',glc.iyear)
    return
end

for imonth = 1:glc.nmonth
    fprintf('Physio-EnergyBal: year=%04d,imonth=%02d\n',glc.iyear,imonth)
    LAI_this_month = glc.LAI_this_site_this_year(imonth);
    if LAI_this_month<0.05 || isnan(LAI_this_month)
        fprintf('Physio-EnergyBal, Skipping month (year=%04d,imonth=%02d), LAI is not avaliable (NaN or below 0.05)\n',glc.iyear,imonth)
        continue % make sure LAI has to be at least 0.05 or cannot be NaN
    end
    
    %% get this month flux data
    [glc,vis_hourly,sw_hourly,lw_hourly,...
    TA_F,SW_IN_F,PA_F,ea_F,qa_F,SWC_F,USTAR_F,ca,LE_F,...
    day_mask,GPP_F]=get_this_month_fluxdata(this_year_data,rad_data,imonth,ca_this_year,glc,eeo);
    
    %% mask nighttime this month flux data at hourly
    TA_F(~day_mask)=nan;
    SW_IN_F(~day_mask)=nan;
    PA_F(~day_mask)=nan;
    SWC_F(~day_mask)=nan;
    LE_F(~day_mask)=nan;
    USTAR_F(~day_mask)=nan;
    qa_F(~day_mask)=nan;
    GPP_F(~day_mask)=nan;
    
    num_invalid_input(1)=sum(isnan(SWC_F));
    num_invalid_input(2)=sum(isnan(TA_F));
    num_invalid_input(3)=sum(isnan(SW_IN_F));
    num_invalid_input(4)=sum(isnan(PA_F));
    num_invalid_input(5)=sum(isnan(qa_F));
    num_invalid_input(6)=sum(isnan(USTAR_F));
    num_invalid_input(7)=sum(isnan(lw_hourly(1,:)));
    
    if num_invalid_input(1)==glc.nHour || num_invalid_input(2)==glc.nHour || num_invalid_input(3)==glc.nHour ||...
            num_invalid_input(4)==glc.nHour || num_invalid_input(5)==glc.nHour || num_invalid_input(6)==glc.nHour...
            || num_invalid_input(7)==glc.nHour || length(unique(num_invalid_input))>1
        fprintf('Physio-EnergyBal, Skipping month (year=%04d,imonth=%02d).  At least one of the followings are missing, or have inconsistent number of time stamps (SWC, Ta, Pz, SW, qa, USTAR, LW_abs)\n',glc.iyear,imonth)
        continue
    end
    
    rb_F = 1./seb.Cv .*(USTAR_F./seb.dleaf).^(-0.5); % CLM5 technical note p60
    rhoa_F = (PA_F - ea_F)./(seb.R.*TA_F) + ea_F./(seb.Rv.*TA_F); % rho air
    
    %% monthly acclimation, need to evaluate Vcmax every month. For other types, this will be skipped.
    if acclimation_type==1
        opt_Vcmax=Vcmax_yearly_peak_month(this_year_data,rad_data,...
            imonth,carbon,seb,glc,eeo,ca_this_year,beta0,within_day_scale);
        opt_Vcmax.leaf_flush_type = leaf_flush_type;
    end
    
    if opt_Vcmax.valid==false
        fprintf('Physio-EnergyBal, Skipping month (year=%04d,imonth=%02d), Opt Vcmax is not avaliable (check peak LAI month meteorological forcing if acclimation_type=2, or check forcing of the current month if acclimation_type=1)\n',glc.iyear,imonth)
        continue
    end
    
    %% start energy balance iteration, hourly
    %opt_Vcmax.hr_daytime_mask = 2; %  hourly
    opt_Vcmax.acclimation_type =acclimation_type; % if 1, acclimate monthly; if 2, acclimate to peak month of each year
    
    valid_mask_hr = ~isnan(vis_hourly);
    Ta_hr = repmat(TA_F',[eeo.Nlayers,1]);Ta_hr(~valid_mask_hr)=nan;
    Pz_hr = repmat(PA_F',[eeo.Nlayers,1]);Pz_hr(~valid_mask_hr)=nan;
    qa_hr = repmat(qa_F',[eeo.Nlayers,1]);qa_hr(~valid_mask_hr)=nan;
    SWC_hr= repmat(SWC_F',[eeo.Nlayers,1]);SWC_hr(~valid_mask_hr)=nan;
    rb_hr = repmat(rb_F',[eeo.Nlayers,1]);rb_hr(~valid_mask_hr)=nan;
    rhoa_hr = repmat(rhoa_F',[eeo.Nlayers,1]);rhoa_hr(~valid_mask_hr)=nan;
    rb_hr_v = rb_hr./(seb.Dv_Dh)^0.67;rb_hr_v(~valid_mask_hr)=nan;
    gb_hr_c = 1./rb_hr_v ./seb.a_CO2_rb;gb_hr_c(~valid_mask_hr)=nan;

    GPP_hr= GPP_F';
    
    G=0;% assume no heat storage in leaves
    albedo_hr=0;% already considered in sw input
    Tleaf_hr = Ta_hr; % for the first guess, assume Tleaf=Ta
    
    max_nvalid_layer = max(sum(~isnan(vis_hourly),1));
    SWC_temp = SWC_hr(1,:).* repmat(linspace(glc.swc_min_pct,glc.swc_max_pct,max_nvalid_layer)',[1,glc.nHour]);
    SWC_hr=nan(eeo.Nlayers,glc.nHour);SWC_hr(1:max_nvalid_layer,:)=SWC_temp;
    
    for itr = 1: eeo.nitr
        [A,g_mol,gs_CO2,Vcmax,Jmax,J,Tleaf_opt,ci,lambda,WUE,Tr,PPFD_chl, D_leafair, Rubisco_idx,ref_vcmax]...
            =Canopy_photosynthesis(Ta_hr,Tleaf_hr,Pz_hr,qa_hr,SWC_hr,ca,vis_hourly/eeo.DeltaL*eeo.cL,gb_hr_c,opt_Vcmax,beta0,carbon,eeo);
              
        output=TRM_second(sw_hourly./eeo.DeltaL, lw_hourly./eeo.DeltaL./seb.emis_l, Ta_hr, qa_hr, albedo_hr, ...
            rb_hr,rb_hr_v,1./(carbon.a.*gs_CO2), seb.emis_l, rhoa_hr, Pz_hr, G, seb);  
        
        delta_Tleaf = output.Ts - Tleaf_hr;
        
        % update Tleaf
        Tleaf_hr = output.Ts;
       
        max_delta_Tleaf=max(abs(delta_Tleaf(:)),[],'omitnan');
        if itr==eeo.nitr
            fprintf('itr=%d (max itr reached),max delta_Tleaf = %f,year=%04d,imonth=%02d\n',itr,max_delta_Tleaf,glc.iyear,imonth);
        end
        if max_delta_Tleaf<eeo.eps_del_Tleaf && itr>1
            break
        end
    end
    
    %% leaf level
    % for global model output
    data2save_hr.leaf.A(:,:,imonth)=A; % μmol CO2 m^-3 s^-1(per unit LAI)
    data2save_hr.leaf.Tr(:,:,imonth)=Tr; % mol H2O m^-3 s^-1 (per unit LAI)
    data2save_hr.leaf.gsCO2_mol(:,:,imonth)=g_mol;% two-side leaf, mol CO2 m^-3 s^-1 (per unit LAI)
    data2save_hr.leaf.ci(:,:,imonth)=ci;% μmol CO2 mol^-1 air
    data2save_hr.leaf.lambda(:,:,imonth)=lambda; % μmol CO2 mol^-1 H2O, (~on the order of 10^3)
    data2save_hr.leaf.PPFD_abs_chl(:,:,imonth) = PPFD_chl; %  mol photon m^-2 s^-1, to get W m^-3 (per unit LAI),divide by 4.55
    
    data2save_hr.leaf.Tleaf4A(:,:,imonth)=Tleaf_opt; % K
    data2save_hr.leaf.H(:,:,imonth) = output.H; % W m^-3 (per unit LAI)
    data2save_hr.leaf.LE(:,:,imonth) = output.LE; % W m^-3 (per unit LAI)
    
    % for internal analysis use
    data2save_hr.leaf.PAR_abs(:,:,imonth) = vis_hourly/eeo.DeltaL; % W m^-3 (per unit LAI) (total)
    data2save_hr.leaf.SW_abs(:,:,imonth) = sw_hourly./eeo.DeltaL; % W m^-3 (per unit LAI), PAR + NIR
    data2save_hr.leaf.LW_abs(:,:,imonth) = lw_hourly./eeo.DeltaL; % W m^-3 (per unit LAI)
   
    data2save_hr.leaf.Vcmax(:,:,imonth)=Vcmax; % at operating temperature
    data2save_hr.leaf.Jmax(:,:,imonth)=Jmax; % at operating temperature
    data2save_hr.leaf.J(:,:,imonth)=J; % at operating temperature
    data2save_hr.leaf.gsCO2(:,:,imonth)=gs_CO2; % two-side leaf
    
    data2save_hr.leaf.Vcmax_m(:,imonth) = ref_vcmax.Vcmax_m; % at reference optimal temperature
    data2save_hr.leaf.Jmax_m(:,imonth) = ref_vcmax.Jmax_m; % at reference optimal temperature
    data2save_hr.leaf.T_m(:,imonth) = ref_vcmax.T_m; % the reference optimal temperature
    if strcmp(within_day_scale,'hourly') % Vcmax acclimate to different hour within a day
        data2save_hr.leaf.hour_idx(:,imonth) = opt_Vcmax.hour_idx; % hour idx
    end  
    data2save_hr.leaf.Rn(:,:,imonth) = output.Rn; % W m^-3 (per unit LAI)
    data2save_hr.leaf.energy_balance(:,:,imonth)=output.energy_balance; % W m^-3 (per unit LAI)
    
    %% canopy level
    % for global model output
    data2save_hr.cano.GPP(:,imonth)=sum(A,1,'omitnan') * eeo.DeltaL; % μmol CO2 m^-2 s^-1
    data2save_hr.cano.Tr(:,imonth)=sum(Tr,1,'omitnan') * eeo.DeltaL; % mol H2O m^-2 s^-1
    data2save_hr.cano.gsCO2_mol(:,imonth)=sum(g_mol,1,'omitnan') * eeo.DeltaL ;% two-side leaf, mol CO2 m^-2 s^-1

    data2save_hr.cano.PPFD_abs_chl(:,imonth) = sum(PPFD_chl,1,'omitnan'); % mol photon m^-2 s^-1, to get W m^-3 (per unit LAI),divide by 4.55
    data2save_hr.cano.FPAR(:,imonth) = rad_data.FAPAR_hourly(:,imonth); % unitless
    
    % for internal analysis use
    data2save_hr.cano.Vcmax(:,imonth)=mean(Vcmax,1,'omitnan'); % μmol CO2 m^-2 s^-1
    data2save_hr.cano.Jmax(:,imonth)=mean(Jmax,1,'omitnan'); % μmol electron m^-2 s^-1
    data2save_hr.cano.J(:,imonth)=mean(J,1,'omitnan'); % μmol electron m^-2 s^-1
    
    data2save_hr.cano.Vcmax_wA(:,imonth)=sum(ci.*Vcmax,1,'omitnan')/sum(A,1,'omitnan'); % μmol CO2 m^-2 s^-1
    data2save_hr.cano.Jmax_wA(:,imonth)=sum(Jmax.*A,1,'omitnan')/sum(A,1,'omitnan'); % μmol electron m^-2 s^-1
    data2save_hr.cano.J_wA(:,imonth)=sum(J.*A,1,'omitnan')/sum(A,1,'omitnan'); % μmol electron m^-2 s^-1
  
    %% input data at canopy level
    data2save_hr.inputs.DeltaL = eeo.DeltaL;
    data2save_hr.inputs.fchl = eeo.cL;

    data2save_hr.inputs.LAI(1,imonth)=LAI_this_month; % from MODIS   
    data2save_hr.inputs.maxLAImonth(1,imonth)=max_month_this_year;
    data2save_hr.inputs.maxLAI(1,imonth)=max_LAI_this_year;
    
    data2save_hr.inputs.SW(:,imonth)=SW_IN_F(1,:);   
    data2save_hr.inputs.Ta(:,imonth)=Ta_hr(1,:);
    data2save_hr.inputs.Pz(:,imonth)=Pz_hr(1,:);
    data2save_hr.inputs.qa(:,imonth)=qa_hr(1,:);
    data2save_hr.inputs.SWC(:,:,imonth)=SWC_hr(:,:);
    data2save_hr.inputs.ca(1,imonth)=ca;
    
    data2save_hr.inputs.GPP_flux(:,imonth)=GPP_hr;

end
% get normalized Vcmax by temp
[output]=normalize_Vcmax(data2save_hr,carbon,eeo); % get normalized Vcmax
data2save_hr.leaf.Vcmax_25 = output.Vcmax25_leaf; % at 25 deg, μmol C m^-3 s^-1 (per unit LAI)
data2save_hr.leaf.Jmax_25 = output.Jmax25_leaf; % at 25 deg, μmol electron m^-3 s^-1 (per unit LAI)
data2save_hr.cano.Vcmax_25_wA = output.Vcmax25_canoA; % μmol CO2 m^-2 s^-1
data2save_hr.cano.Jmax_25_wA = output.Jmax25_canoA; % μmol electron m^-2 s^-1

data2save_hr.cano.Vcmax_25 = output.Vcmax25_cano;
data2save_hr.cano.Jmax_25 = output.Jmax25_cano;

end

