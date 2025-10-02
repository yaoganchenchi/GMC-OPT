function [ qsat ] = qs( T, P )
R   = 287.058;       % dry air gas constant J/kg/K
Rv  = 461.5;       % water vapor gas constant J/kg/K

globalconstant.epsilon= R/Rv;

esat = 611*exp(17.27*(T-273.15)./(T-273.15+237.3));  % Eq. 3.9a from Dingman in Pa

qsat = globalconstant.epsilon*esat./P; % not account for vapor pressure 

end
