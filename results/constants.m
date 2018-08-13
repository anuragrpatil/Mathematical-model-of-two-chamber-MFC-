function const = constants()
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

rmax = 4.2/(60);
ks = 1.27*10^-3;
Yac = 0.212;
bina = 0.02/(24*60);
ds = 0.2544*10^(-4)/(24*60);
rho = 50 ;
dco2 = 0.9960*10^-4/(24*60);
dh = 2.3328*10^-4/(24*60);
Ll = 2*10^(-4);
R = 8.314;
F = 96450; 
T = 303;
Am = 54*10^-4;
bdet = 0.05/(24*60);
kla = 414/(60*24);
co2equi = 7.26*10^-3;
qo2 =2.64/(60*24);
E0cathode = 340*10^(-3);
E0anode =40*10^(-3);
n = 12;
m = 4;
b = 120*10^(-3);
i0 = 0.001*10^(-3);
km =1.7*10^-3;
dm = 4.5;
dcell = 2.5*10^-2;
kaq = 3500*10^(-3);