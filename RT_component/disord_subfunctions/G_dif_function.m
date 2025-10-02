function Gdif=G_dif_function(ng,xg,wg,gL,hL,Gdif)
% this routine evaluates the G function in all quadrature directions
% and checks for normalization

%------Variables ------
%   xg(ng), wg(ng)
% 	gL(ng), hL(ng)
%   phiprime, muprime
% 	upperlimit_pp, lowerlimit_pp, conv1_pp, conv2_pp
% 	Gdif(ng,ng), sum_tp, sum_pp

%  conversion factors to have the ordinates simulate phiprime
upperlimit_pp = 2.0*pi;
lowerlimit_pp = 0.0;
conv1_pp = (upperlimit_pp-lowerlimit_pp)/2.0;
conv2_pp = (upperlimit_pp+lowerlimit_pp)/2.0;

%  now get the Gdif matrix direction by direction
for i = 1: ng
    muprime = xg(i);
    for j = 1: ng
        phiprime = conv1_pp*xg(j) + conv2_pp;
        Gdif(j,i)=G_dir_function(ng,xg,wg,gL,hL,muprime,phiprime);
    end
end

%  check for normalization
%  (1/2PI) int_0^2PI dphi^prime int_0^1 dmu^prime G(OMEGA^prime) = 0.5
%{
sum_tp = 0.0;
for i = (ng/2)+1 : ng
    sum_pp = 0.0;
    for j = 1 : ng
        sum_pp = sum_pp + wg(j)*Gdif(j,i);
    end
    sum_pp = sum_pp*conv1_pp;
    sum_tp = sum_tp + wg(i)*sum_pp;
end
sum_tp = sum_tp/(2.0*pi);
%fprintf('Gfun check (=0.5?): %.8f\n',sum_tp);
%}

end

