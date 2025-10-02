function I_d_uc_d=I_d_uncol_down(DeltaL,Nlayers,Gdif,I_d,ng,xg,I_d_uc_d)
% this routine evaluates the downward uncollided diffuse sky radiation

%------Variables ------
% 	DeltaL, Gdif(ng,ng), I_d, xg(ng)
% 	L1, L2, Prob
%   I_d_uc_d(Nlayers,ng,ng)

% evaluate downward uncollided diffuse sky intensity layer by layer in all
% downward directions
for i = 1: ng/2
    for j = 1: ng
        L1 = 0.0;
        % for k
        k= 1:Nlayers;
        L2 = k*DeltaL - (0.5*DeltaL);
        Prob = exp(- (1/abs(xg(i))) * Gdif(j,i) * (L2-L1) );
        I_d_uc_d(:,j,i) = I_d * Prob;
        % end for k
    end
end
end

