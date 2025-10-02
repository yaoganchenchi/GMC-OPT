function [esat] = es(T)
esat = 611*exp(17.27*(T-273.15)./(T-273.15+237.3));  % Eq. 3.9a from Dingman in Pa
end
