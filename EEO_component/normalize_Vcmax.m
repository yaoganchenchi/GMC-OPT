%% get Vcmax 25
function [output]=normalize_Vcmax(data,carbon,eeo)
% Convert Vcmax_m to Vcmax_25
temp = squeeze(mean(data.leaf.Vcmax_m,1,'omitnan'));
temp2 = squeeze(mean(data.leaf.T_m,1,'omitnan'));

if ndims(data.cano.GPP)==2  % for hourly data
    data.cano.GPP2 = squeeze(mean(data.cano.GPP,1,'omitnan'));
    imag_mask_cano = imag(data.cano.GPP2) ~= 0 | imag(temp) ~= 0 | imag(temp2) ~= 0;  
    GPP_opt = squeeze(sum(data.leaf.A,2,'omitnan'));
else % for non hourly acclimation
    imag_mask_cano = imag(data.cano.GPP) ~= 0 | imag(temp) ~= 0 | imag(temp2) ~= 0;
    GPP_opt = data.leaf.A;
end
imag_mask_leaf =repmat(imag_mask_cano,[1,1,eeo.Nlayers]);
imag_mask_leaf = permute(imag_mask_leaf,[3,1,2]);

% reference Vcmax
GPP_opt(imag_mask_leaf)=nan;
Vcmax_ref = data.leaf.Vcmax_m;Vcmax_ref(imag_mask_leaf)=nan;
Jmax_ref = data.leaf.Jmax_m;Jmax_ref(imag_mask_leaf)=nan;
T_m = data.leaf.T_m;T_m(imag_mask_leaf)=nan;

% Vcmax
del_S = carbon.as-carbon.bs.*(carbon.T_25-carbon.T_ref); % intermedia variable
xx = Arrh(carbon.T_25,carbon.del_H,T_m, carbon.R).*(1+exp((T_m.*del_S-carbon.Hd)./(T_m.*carbon.R)))./(1+exp((carbon.T_25.*del_S-carbon.Hd)./(carbon.T_25.*carbon.R))); % intermedia variable
output.Vcmax25_leaf = Vcmax_ref.*xx;  % leaf
output.Vcmax25_canoA = squeeze(sum(output.Vcmax25_leaf.*GPP_opt,1,'omitnan')./sum(GPP_opt,1,'omitnan')); % canopy weighted by GPP
output.Vcmax25_cano = squeeze(mean(output.Vcmax25_leaf,1,'omitnan')); % arithmetic sum
output.Vcmax25_cano(output.Vcmax25_cano==0)=nan;

% Jmax and J
del_Sj = carbon.asj-carbon.bsj.*(carbon.T_25-carbon.T_ref);% intermedia variable
xxj = Arrh(carbon.T_25,carbon.del_Hj,T_m, carbon.R).*(1+exp((T_m.*del_Sj-carbon.Hd)./(T_m.*carbon.R)))./(1+exp((carbon.T_25.*del_Sj-carbon.Hd)./(carbon.T_25.*carbon.R)));% intermedia variable
output.Jmax25_leaf = Jmax_ref.*xxj;
output.Jmax25_canoA = squeeze(sum(output.Jmax25_leaf.*GPP_opt,1,'omitnan')./sum(GPP_opt,1,'omitnan'));
output.Jmax25_cano = squeeze(mean(output.Jmax25_leaf,1,'omitnan'));
output.Jmax25_cano(output.Jmax25_cano==0)=nan;

end