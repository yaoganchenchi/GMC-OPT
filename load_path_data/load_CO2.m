function ca_NOAA = load_CO2(glc,CO2_path)

CO2 = loadMatData(CO2_path);
ca_temp = [CO2{:,3}];
CO2_syear=1700;
id_year_month = (glc.MODIS_first_year-CO2_syear)*glc.nmonth+1 : (glc.flux_last_year-CO2_syear+1)*glc.nmonth;
ca_temp = ca_temp(id_year_month);
ca_NOAA = nan(1,glc.nyear*glc.nmonth);
idx2 = 1:(glc.flux_last_year-glc.MODIS_first_year+1)*glc.nmonth;
ca_NOAA(idx2) = ca_temp;

end