function Ic=SWEEP_DOWN(Nlayers,ng,xg,wg,Gdif,DeltaL,JJ,R_s,Ic)
% this routine sweeps downwards in the phase-space mesh and handle the
% bottom boundary condition and evaluate the upward Ic at the ground

%------Variables ------

% 	xg(ng), wg(ng)
% 	Gdif(ng,ng), DeltaL, JJ(Nlayers,ng,ng), R_s
% 	fij, aij, bij
% 	Ic(Nlayers+1,ng,ng)
% 	upperlimit, lowerlimit, conv, sum1, sum2
% 	Fc_soil

% sweep downwards
for i = 1: ng/2
    for j = 1: ng
        fij = (xg(i)/DeltaL) - (0.5*Gdif(j,i));
        aij = ((0.5*Gdif(j,i)) + (xg(i)/DeltaL)) / fij;
        bij = 1.0/fij;
        for k = 1: Nlayers
            Ic(k+1,j,i) = aij*Ic(k,j,i) - bij*JJ(k,j,i);
        end
    end
end

% evaluate flux density incident on the ground
sum1 = 0.0;
for i = 1: ng/2
    upperlimit = 2.0*pi;
    lowerlimit = 0.0;
    conv = (upperlimit-lowerlimit)/2.0;
    sum2 = 0.0;
    for j = 1: ng
        sum2 = sum2 + wg(j)*Ic(Nlayers+1,j,i);
    end
    sum1 = sum1 + wg(i)*abs(xg(i))*sum2*conv;
end
Fc_soil = sum1;

% evluate Ic upward at the ground
for i = (ng/2)+1: ng
    for j = 1: ng
        Ic(Nlayers+1,j,i) = (Fc_soil*R_s)/pi;
    end
end

end

