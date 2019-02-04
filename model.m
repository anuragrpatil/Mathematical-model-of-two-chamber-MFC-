% function [ch] = model(T)
clc
clear all 
close all 

T = 273 %:10:453;  
d_t=1;
t = 0:d_t:18000;
cs_temp = zeros(length(T),length(t));

Ds=(0.2544*10^(-4))/86400;  %diffusivity of lactate (m2/day)
Ll=2*10^(-4);               %the thickness of laminar diffusion sublayer()
rho=50;                     %The biomass density (kg/m3)
Yac=0.212;                  %The bacterial yield (kg dry cell/kg)
bina=0.02/86400;            %Inactivation coefficient (/day) 
bdet=0.05/86400;            %Detachment coefficient(/day)
Dco2=(0.9960*10^(-4))/86400; %The diffusivity of carbon dioxide (m2/day)
Dh=(2.3328*10^(-4))/86400;   %Diffusivity of hydrogen(m2/day)
Am=54*10^(-4);               %area of membrane(m2)
F=96450;                    %Fradays constant(C/mol) 
T = T;                    %Temperature(K)
R=8.314;                    %gas constant 
Vc=135*10^(-6);             %
Va=135*10^(-6);             %
Eoanode=340;                % Standard voltage anode (mV)
Eocathode=1299;             % Standard volatage cathode(mV)
ioref=0.001;                % Exchange current in reference conditions(mA)
b=120;                      %Tafel Coefficient (mV)
dcell=2.5*10^(-2);          %distance between electrodes(m)
kaq=3500;                   %Solution conductivity (mS/m)
Eka=-155;                   %
dm=4.5;                     %Membrance thickness(m)
km=1.7;                      %Membrance conductivity (mS/m )
co2equi=7.26*10^(-3);       %
kla=414/86400;              %Overall volumetric oxygen mass tranfer coefficient (/day)
qo2=2.64/86400;             %specific uptake rate of oxygen (/day)
rmax=100.9/3600;            % reaction rate constant(/day)
ks=1.27*10^(-3);            % Half-max-rate lactate concentration 
chb=10^(-7);                %  Concentration of hydrogen in bulk liquid(kg/m3)





f1 = @(cs) rmax*(cs/(ks+cs)) ;
f2 = @(mu,nact,phia) (mu*(1/(1+exp(-F/(R*T*1000)*nact)))*phia) ;
f3 = @(rs,phia,delta,L) Yac*rs-bina*phia+(phia/L)*delta-(phia/L)*(Yac*rs*L+delta) ;
f4 = @(rs,L,delta) Yac*rs*L+delta ;
f5 = @(L) -bdet*L ;
f6 = @(L,csb,cs,rs,delta) (Ds/(Ll*L))*(csb-cs)-rho*rs-(cs/L)*(Yac*rs*L+delta);
f7 = @(L,cco2b,cco2,rs,delta) (Dco2/(Ll*L))*(cco2b-cco2)+4*rho*rs-(cco2/L)*(Yac*rs*L+delta) ;
f8 = @(L,chb,ch,rs,delta) (Dh/(Ll*L))*(chb-ch)+12*rho*rs-(ch/L)*(Yac*rs*L+delta) ;
f9 = @(rs,L,delta) -Am*(Yac*rs*L+delta) ;
f10 = @(csb,cs,vl) (1/vl)*(-(Am*Ds/Ll)*(csb-cs)) ;
f11 = @(cco2b,cco2,vl) (1/vl)*(-(Am*Dco2/Ll)*(cco2b-cco2)) ;
f12 = @(i) -2.24*10^(-7)*(i)+0.0067;
%f12 = @(i) (kla*co2equi+7.7167*10^(-8)*exp(-(i)*(kla+qo2)))/(kla+qo2);
f13 = @(co2,ch)  Eocathode-((R*T*1000)/(4*F))*log(1/(co2*(ch^4))) ;
f14 = @(ch,cco2,cs) Eoanode-((R*T*1000)/(12*F))*log(((cco2^3)*(ch)^(11))/cs) ;
f15 = @(csb) (12*F*Ds*1000*csb)/(Ll) ;
f16 = @(I) ((dm/km)+(dcell/kaq))*(I/1000) ;
f17 = @(I,il) ((R*T)/(12*F))*1000*log(il/(il-I)) ;
f18 = @(I,cs) (b/2.303)*asinh(I/(2*ioref*cs)) ;
f19 = @(Ecathode,Eanode,nohm,nconc,nact) abs(Ecathode-Eanode)-nohm-nconc-nact ;
f20 = @(Eoutput)  Eoutput/100;


d_t=1;
t = 0:d_t:18000;

mu=zeros(1,length(t));
rs=zeros(1,length(t));
phia=zeros(1,length(t));
L=zeros(1,length(t));
cs=zeros(1,length(t));
delta=zeros(1,length(t));
cco2=zeros(1,length(t));
ch=zeros(1,length(t));
vl=zeros(1,length(t));
csb=zeros(1,length(t));
cco2b=zeros(1,length(t));
co2=zeros(1,length(t));
Ecathode=zeros(1,length(t));
Eanode=zeros(1,length(t));
il=zeros(1,length(t));
nohm=zeros(1,length(t));
nconc=zeros(1,length(t));
nact=zeros(1,length(t));
Eoutput=zeros(1,length(t));
I=zeros(1,length(t));

mu(1)=1.166370409*10^(-3); %specific growth rate ()
rs(1)=2.5*10^(-20);        %INITIAL I NEED TO BE KNOWN(TAKEN ZERO HERE) to vary
phia(1)=0.4286802857;      %volume fraction of active biomass 
L(1)=0.025810^(-5);        %Thickness of biofilm 
cs(1)=7.5;                 %The concentration of lactate in biofilm(kg/m3)
delta(1)=-1.446759259*10^(-08); %detachment rate
cco2(1)=1.804718175*10^(-10);%INITIAL CO2 CONC NEEDED and here up needed to be reviewed
ch(1)=10^(-7);               % concentration of hydrogen in biofilm(kg/m3)
vl(1)=135*10^(-6);           %Bulk liquid anoide chamber volume(m3)
csb(1)=5;                    % Concentration of lactate in bulk liquid(kg/m3)
cco2b(1)=1.804718175*10^(-25); %The concentration of hydrogen in bulk of liquid(kg/m3)
co2b(1)=7.23*10^(-3);          % COncnetraiton of oxygen in bulk 
co2(1)=7.23*10^(-3);            %The concentration of oxygen in cathode camber(kg/m3)
Ecathode(1)=845.8563841;        %
Eanode(1)=845.8563841;          %
il(1)=8519.75;                  %limited current of mass transfer 
nohm(1)=0;                     % Ohmic overpontential (mV)
nconc(1)=0;                    %COncentration overpotential(mV)
nact(1)=0;                     %Activation of overpotential(mV)
Eoutput(1)=0;
I(1)=0;

     for i=1:(length(t)-1)

        phia(i+1) = phia(i) + d_t*f3(rs(i),phia(i),delta(i),L(i));
        L(i+1) = L(i) + d_t*f4(rs(i),L(i),delta(i));
        delta(i+1) = delta(i) + d_t*f5(L(i));
        cs(i+1) = cs(i) + 0.00031*d_t*f6(L(i),csb(i),cs(i),rs(i),delta(i));
        cco2(i+1) = cco2(i) + d_t*f7(L(i),cco2b(i),cco2(i),rs(i),delta(i));
        ch(i+1) = ch(i) + 0.0000003*d_t*f8(L(i),chb,ch(i),rs(i),delta(i));
        vl(i+1) = vl(i) + d_t*f9(rs(i),L(i),delta(i));
        csb(i+1) = csb(i) + d_t*f10(csb(i),cs(i),vl(i));
        cco2b(i+1) = cco2b(i) + d_t*f11(cco2b(i),cco2(i),vl(i));
        co2(i+1) = f12(i+d_t);
        Ecathode(i+1) = f13(co2(i+1),ch(i+1));
        Eanode(i+1) = f14(ch(i+1),cco2(i+1),cs(i+1));
        il(i+1) = f15(csb(i+1));
        nohm(i+1) = f16(I(i));
        nconc(i+1) = f17(I(i),il);
        nact(i+1) = f18(I(i),cs(i+1));
        mu(i+1) = f1(cs(i+1));
        rs(i+1) = f2(mu(i+1),nact(i+1),phia(i+1));
        Eoutput(i+1) = f19(Ecathode(i+1),Eanode(i+1),nohm(i+1),nconc(i+1),nact(i+1));
        I(i+1) = f20(Eoutput(i+1));

 
     end
    
idensity=(I/54)*10^(-3);
% mut=(mu(1)-kd).*t;
% mutnoexp(:,:)=exp(mut);
power=Eoutput.*I/10;
eff(:,:)=power/76.7;     
     
% % end
% figure (3)
% % plot(cs(1:14000),eff(1:14000))
% plot(eff(1:12000)/10,cs(1:12000))
% xlabel('Efficiency')
% ylabel('Concentration of lactate (kg/m^{3})')
% saveas(gcf,'efficiency.tiff')
% % xlim([2.1,7.5]);
% % ylim([]);
% figure (4)
% plot(I(1:14000)/10,cs(1:14000))
% xlabel('Current density (mA/m^{2})')
% ylabel('Concentration of lactate (kg/m^{3})')
% saveas(gcf,'current_density.tiff')
% % xlim([]);
% % ylim([]);



% figure (7)
% plot(Eanode(11:14000),I(11:14000),'k','LineWidth',0.5)
% xlabel('anode potential(mV)','FontWeight','bold')
% ylabel('Current(mA)','FontWeight','bold')
% saveas(gcf,'current_AnodePotential.tiff')

figure (8)
plot(Eoutput(1:14000)/10,idensity(1:14000)/10,'k','LineWidth',1)
xlabel('Voltage(mV)','FontWeight','bold')
ylabel('Current denisty(mA/mm{2})','FontWeight','bold')
% xlim([15,35]);
% ylim([2.5,7.55]);
saveas(gcf,'currentdensityVSvoltage.tiff')


figure (9)
plot(cs(11:14000),idensity(11:14000),'k','LineWidth',1)
xlabel('Concentration of lactate (kg/m^{3})','FontWeight','bold')
ylabel('Current denisty(mA/mm^{2})','FontWeight','bold')
xlim([3,7.3]);
% ylim([2.5,7.55]);
saveas(gcf,'currentdensityVScs.tiff')



