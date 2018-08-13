function [Ddv_Div]=DEdef(I,D)

Ds=(0.2544*10^(-4))/86400;
Ll=2*10^(-4);
rho=50;
Yac=0.212;
bina=0.02/86400;
bdet=0.05/86400;
Dco2=(0.9960*10^(-4))/86400;
Dh=(2.3328*10^(-4))/86400;
Am=54*10^(-4);
F=96450;
T=303;
R=8.314;
Vc=135*10^(-6);
Va=135*10^(-6);
Eoanode=340;
Eocathode=1299;
ioref=0.001;
b=120;
dcell=2.5*10^(-2);
kaq=3500;
Eka=-155;
dm=4.5;
km=1.7;
co2equi=7.26*10^(-3);
kla=414/86400;
qo2=2.64/86400;
rmax=100.9/3600;
ks=1.27*10^(-3);
chb=10^(-7);

mu=D(1);
rs=D(2);
phia=D(3);
L=D(4);
cs=D(5);
delta=D(6);
cco2=D(7);
ch=D(8);
vl=D(9);
csb=D(10);
cco2b=D(11);
co2=D(12);
Ecathode=D(13);
Eanode=D(14);
il=D(15);
nohm=D(16);
nconc=D(17);
nact=D(18);
Eoutput=D(19);
i=D(20);


Ddv_Div=[ rmax*(cs/(cs+ks));%1
          mu*(1/(1+exp((-F/(R*T))*(nact/1000))))*phia;%2
          Yac*rs-bina*phia+(phia/L)*delta-(phia/L)*(Yac*rs*L+delta);%3
          Yac*rs*L+delta;%4
          (Ds/(Ll*L))*(csb-cs)-rho*rs-(cs/L)*(Yac*rs*L+delta);%5
          -bdet*L;%6
          (Dco2/(Ll*L))*(cco2b-cco2)+4*rho*rs-(cco2/L)*(Yac*rs*L+delta);%7
          (Dh/(Ll*L))*(chb-ch)+12*rho*rs-(ch/L)*(Yac*rs*L+delta);%8
          -Am*(Yac*rs*L+delta);%9
          (1/vl)*((-Am*Ds/Ll)*(csb-cs));%10
          (1/vl)*((-Am*Dco2/Ll)*(cco2b-cco2));%11
          kla*(co2equi-co2)-qo2*co2;%12
          Eocathode-((R*T*1000)/(4*F))*log(1/(co2*(ch^4)));%13
          Eoanode-((R*T*1000)/(12*F))*log(((cco2^3)*(ch)^(11))/cs);%14
          (12*F*Ds*1000*csb)/(Ll);%15
          ((dm/km)+(dcell/kaq))*(i/1000);%16
          ((R*T)/(12*F))*1000*log(il/(il-i));%17
          (b/2.303)*asinh(i/(2*ioref*cs));%18
          abs(Ecathode-Eanode)-nohm-nconc-nact;%19
          Eoutput/100];%20
          
          
       
end
          
          
          
    

