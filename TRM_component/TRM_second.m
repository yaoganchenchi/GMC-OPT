function [output] = TRM_second(SWin, LWin, Ta, q, albedo, rah,raw, rs, emis, rhoa, Ps, G, globalconstant)

%% deal with ground heat flux and its derivatives with respect to albedo, ra and rs
Rn_star = SWin.*(1-albedo) + emis.*LWin - emis.*globalconstant.sb.*(Ta).^4; % potential to acutual temperature

%% important intermediate variables
delta   = desdT(Ta);
gamma = (globalconstant.cp.*Ps)./(globalconstant.epsilon.*globalconstant.Lv); %(globalconstant.cp.*Ps)./(0.622.*globalconstant.Lv);
qsat    = qs(Ta, Ps);
lambda_o_prime = 1./(rhoa*globalconstant.cp);
ro =  rhoa*globalconstant.cp./(4*emis.*globalconstant.sb.*Ta.^3);
f_TRM = 1./ro+(1./rah).*(1 + (delta./gamma).*(rah./(raw + rs)));
AA = rhoa.*globalconstant.Lv.*(qsat-q);

%% computation of Ts, H, LE, and energy balance based on second-order Taylor expansion SEB (this is from Paw U 1987, J. Therm. Biol) 
a = lambda_o_prime.*6.*emis.*globalconstant.sb.*Ta.^2 + 1/2*(des2dT2(Ta)./gamma).*(1./(raw+rs)); 
b = f_TRM;
c = -lambda_o_prime.*(Rn_star - G -AA./(raw+rs));

Ts1 = Ta + (-b + sqrt(b.^2 - 4*a.*c))./(2*a);

Ts = Ts1;
H = rhoa.*globalconstant.cp.*(Ts - Ta)./rah;
LE = rhoa.*globalconstant.Lv.*(qs(Ts, Ps)-q)./(raw+rs); 
energy_balance = SWin.*(1-albedo) + emis.*LWin - emis.*globalconstant.sb.*Ts.^4  - H - LE - G;
Rn = SWin.*(1-albedo) + emis.*LWin - emis.*globalconstant.sb.*Ts.^4;

output.Ts=Ts;
output.H=H;
output.LE=LE;
output.energy_balance=energy_balance;
output.Rn=Rn;

end







