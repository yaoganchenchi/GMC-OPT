%% Arrhenius type function
function farrh=Arrh(T,del_H,T_ref, R)
% most T_ref = 298.15K
% R = 8.3145 J mol^-1 K-1
farrh = exp(del_H.*(T-T_ref)./(T_ref.*R.*T));
end
