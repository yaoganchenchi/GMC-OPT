function [I_o_uc_u]=I_o_uncol_up(DeltaL,Nlayers,Gdir,Gdif,I_o,mu_o,ng,xg,R_s,I_o_uc_u)
% this routine evaluates the upward uncollided direct solar radiation

%------Variables ------
% 	DeltaL, Gdir, Gdif(ng,ng), I_o
% 	mu_o, R_s, xg(ng), L1
% 	L2, Prob, I_o_uc_u_soil, F_o_uc_d_soil
%   I_o_uc_u(Nlayers,ng,ng)

% the uncollided direct solar radiation incident on the ground
L1 = 0.0;
L2 = Nlayers*DeltaL;
F_o_uc_d_soil = abs(mu_o)*I_o*exp(-(1/abs(mu_o))*Gdir*(L2-L1));

% upward uncollided intensity from reflection by soil of I_o_uc_d_soil
I_o_uc_u_soil = (R_s/pi)*F_o_uc_d_soil;

% evaluate upward uncollided direct solar intensity layer by layer in all upward directions
for i = (ng/2)+1: ng % has potential to further optimize
    for j = 1: ng
        L2 = Nlayers*DeltaL;
        %for k
        k= 1:Nlayers;
        L1 = k*DeltaL - (0.5*DeltaL);
        Prob = exp(- (1/abs(xg(i))) * Gdif(j,i) * (L2-L1) );
        I_o_uc_u(:,j,i) = I_o_uc_u_soil * Prob;
        %end for k
    end
end
end

