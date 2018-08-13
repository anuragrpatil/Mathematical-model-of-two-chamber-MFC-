clear all
domain=[0 60];

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



IC1=1.166370409*10^(-3);
IC2=2.5*10^(-4);%INITIAL I NEED TO BE KNOWN(TAKEN ZERO HERE) to vary
IC3=0.4286802857;%to vary
IC4=0.025;%to vary
IC5=5;
IC6=-1.446759259*10^(-08);
IC7=1.804718175*10^(-8);%INITIAL CO2 CONC NEEDED and here up needed to be reviewed
IC8=10^(-7);
IC9=135*10^(-6);
IC10=5;
IC11=1.804718175*10^(-8);
IC12=7.23*10^(-3);
IC13=845.8563841;
IC14=845.8563841;
IC15=8519.75;
IC16=0;
IC17=0;
IC18=0;
IC19=0;
IC20=0;


IC=[IC1 IC2 IC3 IC4 IC5 IC6 IC7 IC8 IC9 IC10 IC11 IC12 IC13 IC14 IC15 IC16 IC17 IC18 IC19 IC20];

[IVSOL,DVSOL]=ode23('DEdef',domain,IC);
plot(IVSOL,DVSOL(:,5),'r-')
xlabel('time in hr');
ylabel('cs');

