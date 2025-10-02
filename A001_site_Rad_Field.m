function [SW_rad_field,LW_rad_field] = A001_site_Rad_Field(leaf_angle_distr,vis_albedo_type,eeo,glc,seb,RT_input_this_year_data)
% shortwave and longwave

if ~exist('leaf_angle_distr','var')
    leaf_angle_distr = 'Spherical';
end
if ~exist('vis_albedo_type','var')
    vis_albedo_type = 'vis_median';
end

%% load paths and add dirctory for model tools
A000_load_paths

ng = 6; %Gaussian quadrature, 6 for fast calculation. 12 for more accurate calculation.
rad2deg = 180.0/pi;
deg2rad = pi/180.0;
fdir_lw = 0; % assume all lw radiation is diffuse
SZA_lw = 45; % place holder. Any value is fine. It won't be used in lw rad.

% Initilize
AB_vis_I_stack = nan(eeo.Nlayers,glc.nHour,glc.nmonth);% shortwave
AB_nir_I_stack = nan(eeo.Nlayers,glc.nHour,glc.nmonth);% shortwave

FAPAR_hourly = nan(glc.nHour,glc.nmonth);% shortwave

AB_DOWN_I_stack = nan(eeo.Nlayers,glc.nHour,glc.nmonth);% longwave
AB_UP_I_stack = nan(eeo.Nlayers,glc.nHour,glc.nmonth);% longwave
AB_leaf_source_I_stack = nan(eeo.Nlayers,glc.nHour,glc.nmonth);% longwave

% get number of days in a month
if isLeapYear(glc.iyear)
    day_month = [31,29,31,30,31,30,31,31,30,31,30,31];
else
    day_month = [31,28,31,30,31,30,31,31,30,31,30,31];
end

for imonth = 1:glc.nmonth
    LAI_this_month = glc.LAI_this_site_this_year(imonth);
    if LAI_this_month<0.05 || isnan(LAI_this_month)
        fprintf('Skipping RT: year=%04d,imonth=%02d, LAI is not avaliable (NaN or below 0.05)\n',glc.iyear,imonth)   
        continue % make sure LAI has to be at least 0.05 or cannot be NaN
    else
        fprintf('Running RT: year=%04d,imonth=%02d\n',glc.iyear,imonth)     
    end
    
    nmax_day = day_month(imonth);
    idoy_idx = sum(day_month(1:imonth-1))+1:sum(day_month(1:imonth-1))+nmax_day ;
        
    idoy_count=0;
    sunzenith=nan(glc.nHour,nmax_day);
    sunazimuth=nan(glc.nHour,nmax_day);
    for idoy = idoy_idx
        idoy_count = idoy_count+1;
        for ihour = 1:glc.nHour
            this_hour = ihour*glc.hh_intv-glc.hh_intv;% starting from 00, can consider the middle of the time step: hh_intv/2;
            [sunzenith(ihour,idoy_count),sunazimuth(ihour,idoy_count)]=solor_position_calculator(glc.iyear,idoy,glc.LON,glc.LAT,this_hour);
        end
    end
    SZA_ave = mean(sunzenith,2,'omitnan');
    
    % get this month flux data
    idx_month = (imonth-1)*glc.nrow_per_month+1:imonth*glc.nrow_per_month;
    this_month_data=[];
    fnames = fieldnames(RT_input_this_year_data);
    for ivar = 1:length(fnames)
        this_var = fnames{ivar};
        this_month_data.(this_var) = RT_input_this_year_data.(this_var)(idx_month,:);
    end
    
    % shortwave
    SW_IN_F = this_month_data.SW_IN_F; % shortwave
    SW_IN_POT = this_month_data.SW_IN_POT;% shortwave
    PA_F = this_month_data.PA_F; %kPa
    num_night_day = this_month_data.NIGHT;
     
    LW_IN_F = this_month_data.LW_IN_F; % longwave
    LW_UP_F = (this_month_data.TS_F_MDS_1+273.15).^4.*seb.emis_s.*seb.sb; % longwave
    LW_leaf_source_F  = seb.emis_l.*seb.sb.*(this_month_data.TA_F+273.15).^4 ; % the leaf LW source for each layer

    day_mask = SW_IN_F>10 & num_night_day<10 & SZA_ave < 90; % SZA<90
    
    % get fdir and fdif partioning
    SZA_ave(~day_mask)=nan;
    SZA_ave=  SZA_ave*deg2rad; % go to rad
    SW_IN_F(~day_mask)=nan;  SW_IN_POT(~day_mask)=nan;  PA_F(~day_mask)=nan; % shortwave
    [fdir_vis,fdir_nir]=dir_diffuse_partitioning(SW_IN_F,SW_IN_POT,PA_F,SZA_ave);
   
    SW_IN_VIS_correct = 0.45*SW_IN_F;% shortwave
    SW_IN_NIR_correct = 0.55*SW_IN_F;% shortwave
    SW_IN_VIS = repmat(SW_IN_VIS_correct',[eeo.Nlayers,1]);% shortwave
    SW_IN_NIR = repmat(SW_IN_NIR_correct',[eeo.Nlayers,1]);% shortwave
    
    LW_IN_F(~day_mask)=nan;  LW_UP_F(~day_mask)=nan;     % longwave
    LW_IN = repmat(LW_IN_F',[eeo.Nlayers,1]); % longwave
    LW_UP = repmat(LW_UP_F',[eeo.Nlayers,1]); % longwave
    LW_leaf_source = repmat(LW_leaf_source_F',[eeo.Nlayers,1]); 
    
    % get each hour abs profile, run RT model
    SZA_ave=SZA_ave*rad2deg; % back to deg
    day_hour_idx = find(day_mask'==1);
    AB_vis = nan(eeo.Nlayers,glc.nHour);% shortwave
    AB_nir = nan(eeo.Nlayers,glc.nHour);% shortwave
    
    AB_DOWN = nan(eeo.Nlayers,glc.nHour);% longwave
    AB_UP = nan(eeo.Nlayers,glc.nHour);% longwave
    AB_leaf_source = nan(eeo.Nlayers,glc.nHour);% longwave
    for ihour = day_hour_idx
        % shortwave
        AB=disord1d(180-SZA_ave(ihour),fdir_vis(ihour),LAI_this_month,leaf_angle_distr,vis_albedo_type,ng); % to match the polar coordinate in disord1d
        AB_vis(1:length(AB),ihour)=AB;
        
        AB=disord1d(180-SZA_ave(ihour),fdir_nir(ihour),LAI_this_month,leaf_angle_distr,'nir',ng);
        AB_nir(1:length(AB),ihour)=AB;
        
        % longwave
        AB1=disord1d_lw(180-SZA_lw,fdir_lw,LAI_this_month,leaf_angle_distr,'lw_down',ng);
        AB_DOWN(1:length(AB),ihour)=AB1;
        
        %AB2=disord1d_lw(180-SZA_lw,fdir_lw,LAI_this_month,leaf_angle_distr,'lw_up',ng);
        AB_UP(1:length(AB),ihour)=AB1; % a symmetric problem if assume soil lw abs is 1 
        
        AB2=disord1d_lw_leaf_source(180-SZA_lw,fdir_lw,LAI_this_month,leaf_angle_distr,'leaf_source',ng);
        AB_leaf_source(1:length(AB),ihour)=AB2;
    end
    
    %% shortwave
    AB_vis_I = AB_vis.*SW_IN_VIS;
    AB_nir_I = AB_nir.*SW_IN_NIR;
    AB_vis_I_stack(:,:,imonth) = AB_vis_I;
    AB_nir_I_stack(:,:,imonth) = AB_nir_I;
    
    % FAPAR
    FAPAR_hourly(:,imonth) = sum(AB_vis_I,1,'omitnan')./ SW_IN_VIS_correct';
    
    %% longwave
    n_valid_layer = round(LAI_this_month/eeo.DeltaL);
    temp = flipud(AB_UP(1:n_valid_layer,:));
    AB_UP=nan(size(AB_UP));
    AB_UP(1:n_valid_layer,:) = temp;
    
    AB_DOWN_I = AB_DOWN.*LW_IN;
    AB_UP_I = AB_UP.*LW_UP;
    AB_leaf_source_I = AB_leaf_source.* LW_leaf_source; % assume Telaf = Ta, can also consider no leaf but only air, either way
    AB_DOWN_I_stack(:,:,imonth) = AB_DOWN_I;
    AB_UP_I_stack(:,:,imonth) = AB_UP_I;
    AB_leaf_source_I_stack(:,:,imonth) = AB_leaf_source_I;
    
end

% shortwave
SW_rad_field.AB_vis_I_stack=AB_vis_I_stack;
SW_rad_field.AB_nir_I_stack=AB_nir_I_stack;

SW_rad_field.FAPAR_hourly = FAPAR_hourly;

% longwave
LW_rad_field.AB_lw_stack = AB_DOWN_I_stack + AB_UP_I_stack + AB_leaf_source_I_stack;

end
