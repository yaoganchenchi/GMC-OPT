function [Ic,convergence] = SWEEP_UP(Nlayers,ng,xg,Gdif,DeltaL,JJ,Ic,epsilon)
% this routine sweeps upwards in the phase-space mesh and checks for
% convergence

%------Variables ------
% 	xg(ng), wg(ng), Gdif(ng,ng), DeltaL, JJ(Nlayers,ng,ng)
% 	fij, cij, dij
% 	Ic(Nlayers+1,ng,ng), epsilon, Ic_old(ng,ng)

% save earlier iterate of Ic(1,j,i) for convergence checking
Ic_old = zeros(ng,ng);
Ic_old(1:ng,(ng/2)+1: ng) = Ic(1,1:ng,(ng/2)+1: ng);

% sweep upwards
for i = (ng/2)+1: ng
    for j = 1: ng
        fij = ((xg(i)/DeltaL)+(0.5*Gdif(j,i)));
        cij = ((xg(i)/DeltaL)-(0.5*Gdif(j,i))) / fij;
        dij = 1.0/fij;
        for k = Nlayers:-1: 1
            Ic(k,j,i) = cij*Ic(k+1,j,i) + dij*JJ(k,j,i);
        end
    end
end

% check for convergence
convergence = true;
Ic_dif = abs(squeeze(Ic(1,1:ng,(ng/2)+1: ng)) - Ic_old(1:ng,(ng/2)+1: ng));
if max(Ic_dif(:))> epsilon
    convergence = false;
end

end