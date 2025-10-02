%% load some global constants. Do not change
[carbon,seb,glc,eeo]= getConstants();

carbon.n1=1; % all set to C3
carbon.n2=1; % all set to C3

%%%%%%%%%% Not in use
glc.swc_min_pct=1.;
glc.swc_max_pct=1.;
%%%%%%%%%%

%% specify accoarding to your own data
glc.MODIS_first_year = 2001;
glc.flux_last_year =2014;
glc.nyear = glc.flux_last_year-glc.MODIS_first_year+1;
