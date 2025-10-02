function [I_d_uc_u]=I_d_uncol_up(DeltaL,Nlayers,Gdif,I_d,ng,xg,wg,R_s,I_d_uc_u)
% this routine evaluates the upward uncollided diffuse sky radiation

%------Variables ------
% 	DeltaL, Gdif(ng,ng), I_d, xg(ng), wg(ng), R_s
% 	upperlimit, lowerlimit, conv, sum1, sum2
% 	L1, L2, Prob
%   I_d_uc_d_soil, I_d_uc_u_soil, F_d_uc_d_soil
% 	I_d_uc_u(Nlayers,ng,ng)

% evaluate diffuse sky flux density incident on the soil below canopy
L1 = 0.0;
L2 = Nlayers*DeltaL;
sum1 = 0.0;
for i = 1: ng/2
    upperlimit = 2.0*pi;
    lowerlimit = 0.0;
    conv = (upperlimit-lowerlimit)/2.0;
    sum2 = 0.0;
    for j = 1: ng
        Prob = exp(- (1/abs(xg(i))) * Gdif(j,i) * (L2-L1) );
        I_d_uc_d_soil = I_d * Prob;
        sum2 = sum2 + wg(j)*I_d_uc_d_soil;
    end
    sum1 = sum1 + wg(i)*abs(xg(i))*sum2*conv;
end
F_d_uc_d_soil = sum1;

% evaluate upward uncollided intensity due to reflection from soil
I_d_uc_u_soil = F_d_uc_d_soil*(R_s/pi);

% evaluate downward uncollided diffuse sky intensity layer by layer in all
% upward directions

for i = (ng/2)+1: ng
    for j = 1: ng
        L2 = Nlayers*DeltaL;
        % for k
        k= 1:Nlayers;
        L1 = k*DeltaL - (0.5*DeltaL);
        Prob = exp(- (1/abs(xg(i))) * Gdif(j,i) * (L2-L1) );
        I_d_uc_u(:,j,i) = I_d_uc_u_soil * Prob;
        % end for k
    end
end

end

