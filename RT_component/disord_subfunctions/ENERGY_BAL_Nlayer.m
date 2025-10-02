function [AB_uc,AB_c,AB] = ENERGY_BAL_Nlayer(Nlayers,ng,wg,DeltaL,rho_Ld,tau_Ld,Gdir,Gdif,I_o_uc_d,I_d_uc_d,I_o_uc_u,I_d_uc_u,Ic)
%------Variables ------
% 	xg(ng), wg(ng), mu_o
% 	rho_Ld, tau_Ld, Gdif(ng,ng), DeltaL, R_s
% 	L1, L2, Prob, I_uc_u
% 	I_o_uc_d(Nlayers), I_d_uc_d(Nlayers,ng,ng)
% 	I_o_uc_u(Nlayers,ng,ng), I_d_uc_u(Nlayers,ng,ng)
%  	Ic(Nlayers+1,ng,ng)
% 	upperlimit, lowerlimit, conv, sum1, sum2
% 	AB_c,AB_c,AB_uc

upperlimit = 2.0*pi;
lowerlimit = 0.0;
conv = (upperlimit-lowerlimit)/2.0;

% evaluate canopy absorption from I_o_uc_d
temp1 = I_o_uc_d*Gdir;
AB_o_uc_d = temp1*(1.0-(rho_Ld+tau_Ld))*DeltaL;

% evaluate canopy absorption from I_o_uc_u
AB_o_uc_u=zeros(Nlayers,1);
for k = 1: Nlayers
    sum1 = 0.0;
    for i = (ng/2)+1: ng
        sum2 = 0.0;
        for j = 1: ng
            sum2 = sum2 + wg(j)*I_o_uc_u(k,j,i)*Gdif(j,i);
        end
        sum1 = sum1 + wg(i)*sum2*conv;
    end
    AB_o_uc_u(k) = sum1*(1.0-(rho_Ld+tau_Ld))*DeltaL;
end

% evaluate canopy absorption from I_d_uc_d
AB_d_uc_d=zeros(Nlayers,1);
for k = 1: Nlayers
    sum1 = 0.0;
    for i = 1: ng/2
        sum2 = 0.0;
        for j = 1: ng
            sum2 = sum2 + wg(j)*I_d_uc_d(k,j,i)*Gdif(j,i);
        end
        sum1 = sum1 + wg(i)*sum2*conv;
    end
    AB_d_uc_d(k) = sum1*(1.0-(rho_Ld+tau_Ld))*DeltaL;
end

% evaluate canopy absorption from I_d_uc_u
AB_d_uc_u=zeros(Nlayers,1);
for k = 1:Nlayers
    sum1 = 0.0;
    for i = (ng/2)+1: ng
        sum2 = 0.0;
        for j = 1: ng
            sum2 = sum2 + wg(j)*I_d_uc_u(k,j,i)*Gdif(j,i);
        end
        sum1 = sum1 + wg(i)*sum2*conv;
    end
    AB_d_uc_u(k) = sum1*(1.0-(rho_Ld+tau_Ld))*DeltaL;
end


% evaluate canopy absorption from I_c
AB_c=zeros(Nlayers,1);
for k = 1: Nlayers
    sum1 = 0.0;
    for i = 1: ng
        sum2 = 0.0;
        for j = 1: ng
            sum2 = sum2 + wg(j)*Ic(k,j,i)*Gdif(j,i);
        end
        sum1 = sum1 + wg(i)*sum2*conv;
    end
    AB_c(k) = sum1*(1.0-(rho_Ld+tau_Ld))*DeltaL;
end

AB_uc = AB_o_uc_d + AB_o_uc_u + AB_d_uc_d + AB_d_uc_u;
% AB_up = AB_o_uc_u + AB_d_uc_u;
% AB_down = AB_o_uc_d + AB_d_uc_d;
AB = AB_c + AB_uc;
end