function [S]=MULTI_COLL_S(Nlayers,ng,wg,Gamma_d_dif,Ic,S)
% this routine evaluates the multiple-collision source S

%------Variables ------
% 	xg(ng), wg(ng)
% 	Gamma_d_dif(ng,ng,ng,ng)
% 	upperlimit1, lowerlimit1, conv11, sum1
%	upperlimit2, lowerlimit2, conv21, sum2
% 	Ic(Nlayers+1,ng,ng), Ic_cell_center
% 	S(Nlayers,ng,ng)

%evaluate multiple-collision source due to Ic
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
                    Ic_cell_center = 0.5*(Ic(k,m,n)+Ic(k+1,m,n));
                    sum2 = sum2 + wg(m)*(1.0/pi)*Gamma_d_dif(m,n,j,i)*Ic_cell_center;
                end
                sum1 = sum1 + wg(n)*sum2*conv21;
            end
            S(k,j,i) =  sum1*conv11;
        end
    end
end

end
