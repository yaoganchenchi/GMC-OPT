function I_o_uc_d=I_o_uncol_down(DeltaL,Nlayers,Gdir,I_o,mu_o)
% this routine evaluates the dpwnward uncollided direct solar radiation

%------Variables ------
% 	DeltaL, Gdir, I_o, mu_o
% 	L1, L2, Prob, I_o_uc_d(Nlayers)

% evaluate downward uncollided direct solar intensity layer by layer
L1 = 0.0;
k = 1:Nlayers;
L2 = k*DeltaL - (0.5*DeltaL);
Prob = exp(- (1/abs(mu_o)) * Gdir * (L2-L1) );
I_o_uc_d = I_o * Prob';
end
