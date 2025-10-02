function Q=FCS(Nlayers,ng,wg,Gamma_d_dir,Gamma_d_dif,I_o_uc_d, I_o_uc_u,I_d_uc_d, I_d_uc_u)
% this routine evaluates the first-collision source Q

%------Variables ------
% 	xg(ng), wg(ng)
% 	Gamma_d_dir(ng,ng),Gamma_d_dif(ng,ng,ng,ng)
% 	I_o_uc_d(Nlayers), I_o_uc_u(Nlayers,ng,ng)
% 	I_d_uc_d(Nlayers,ng,ng), I_d_uc_u(Nlayers,ng,ng)
% 	upperlimit1, lowerlimit1, conv11, sum1
% 	upperlimit2, lowerlimit2, conv21, sum2
%  	Q(Nlayers,ng,ng)

% evaluate first-collision source due to I_o_uc_d
I_o_uc_d_mtrx = repmat(I_o_uc_d,[1,ng,ng]);
Gamma_d_dir_mtrx = repmat(Gamma_d_dir,[1,1,Nlayers]);
Gamma_d_dir_mtrx = permute(Gamma_d_dir_mtrx,[3,1,2]);
Q = (1.0/pi).*Gamma_d_dir_mtrx.*I_o_uc_d_mtrx;

% evaluate first-collision source due to I_o_uc_u, I_d_uc_d, I_d_uc_u
for i = 1: ng
    for j = 1: ng
        for k = 1: Nlayers
            upperlimit1 =  1.0;
            lowerlimit1 = -1.0;
            conv11 = (upperlimit1-lowerlimit1)/2.0;
            sum1 = 0.0;
            for n = 1: ng
                upperlimit2 = 2.0*pi;
                lowerlimit2 = 0.0;
                conv21 = (upperlimit2-lowerlimit2)/2.0;
                sum2 = 0.0;
                for m = 1: ng
                    sum2 = sum2 + wg(m)*(1.0/pi)*Gamma_d_dif(m,n,j,i)*(I_o_uc_u(k,m,n)+I_d_uc_d(k,m,n)+I_d_uc_u(k,m,n));
                end
                sum1 = sum1 + wg(n)*sum2*conv21;
            end
            Q(k,j,i) = Q(k,j,i) + sum1*conv11;
            
        end
    end
end

end