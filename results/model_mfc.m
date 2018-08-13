function dy = model_mfc(t,y)
dy = zeros(length(y),1);
global rmax
global ks
global Yac
global bina
global ds
global rho 
global dco2
global dh
global Ll
global R
global F
global T
global Am
global bdet
global kla
global co2equi
global qo2
global E0cathode
global E0anode
global n 
global m 
global b 
global i0
global km 
global dm 
global dcell
global kaq

rs=y(1);
nu=y(2);
phi_a=y(3);
phi_i=y(4);
L=y(5);
delta=y(6);
cs = y(7);
cco2 = y(8);
ch = y(9);
vl = y(10);
csb = y(11);
cco2b = y(12);
co2 = y(13);
Eoutput = y(14);
Ecathode = y(15);
Eanode = y(16);
nact=y(17);
nohm=y(18);
ncon=y(19);
il=y(20);
i=y(21);
chb=y(22);


dy = [
(nu*(1/(1+exp(-F/(R*T)*nact)))*phi_a); %1 al
rmax*(cs/(ks+cs));%2 al
Yac*rs-bina*phi_a+(phi_a/L)*delta-(phi_a/L)*(Yac*rs*L+delta) %3 diff
bina*phi_a+(phi_i/L)*delta-(phi_a/L)*(Yac*rs*L+delta) %4 diff 
Yac*rs*L+delta; %5 diff
-bdet*L; % 6 al
ds*(Ll*L)*(csb-cs)-rho*rs-(cs/L)*(Yac*rs*L+delta); %7 diff
dco2*(Ll*L)*(cco2b-cco2)+4*rho*rs-(cco2/L)*(Yac*rs*L+delta); %8 diff
dh*(Ll*L)*(chb-ch)+12*rho*rs-(cs/L)*(Yac*rs*L+delta); % 9 diff
-Am*L; %10 diff
(1/vl)*(-(Am*ds/Ll)*(csb-cs)); %11 diff
(1/vl)*(-(Am*dco2/Ll)*(cco2b-cco2));% 12 diff
kla*(co2equi-co2)-qo2*co2; % 13 diff
abs(Ecathode-Eanode)-nohm-nact-ncon; %14 al
E0cathode-(R*T/(m*F*2.303))*log(1/(co2*ch^4)); %15 al
E0anode-(R*T/(n*F*2.303))*log(cco2^3*ch^11/cs);%16 al
b/2.303*asinh(i/(2*i0*cs));%17 al
(dm/km+dcell/kaq)*i;%18 al
(R*T/(n*F*2.303))*log(il/(il-i));%19 al
(n*F*ds*csb)/L;%20 al 
Eoutput/100;%21 al
(1/vl)*(-(Am*dh/Ll)*(chb-ch));%22 diff
];



