clear all 
close all 
clc

 
%%%%%%%%%%%%%%%%%%%%%%%%% parameters in the sys of ode%%%%%%%%%%

% Ds=(0.2544*10^(-4))/86400;
Ll=2*10^(-4);
rho=50;
Yac=0.212;
bina=0.02/86400;
bdet=0.05/86400;
% Dco2=(0.9960*10^(-4))/86400;
% Dh=(2.3328*10^(-4))/86400;
Am=54*10^(-4);
F=96450;
T=323;
R=8.314;
Vc=135*10^(-6);
Va=135*10^(-6);
Eoanode=340;
Eocathode=1299;
ioref=0.005;%changed from tafel plot to 0.005 from 0.001
b=120;
dcell=2.5*10^(-2);
kaq=3500;
% Eka=-155;DsEff(i)
dm=4.5;
km=1.7;
co2equi=7.26*10^(-3);
kla=414/86400;
qo2=2.64/86400;
rmax=100.9/3600;
ks=1.27*10^(-3);
chb=10^(-7);


% mL = [0.04*10^(-1), 0.05*10^(-1), 0.065*10^(-1)];
 mT = [321 324 328]

for j=1:(length(mT))
    T = mT(j)
f1 = @(cs) rmax*(1-exp(-cs/ks)) ;
f2 = @(mu,nact,phia) (mu*(1/(1+exp(-F/(R*T*1000)*nact)))*phia) ;
f3 = @(rs,phia,delta,L) Yac*rs-bina*phia+(phia/L)*delta-(phia/L)*(Yac*rs*L+delta) ;
f4 = @(rs,L,delta) Yac*rs*L+delta ;
f5 = @(L) -bdet*L ;
f6 = @(L,csb,cs,rs,delta,Ds) (Ds/(Ll*L))*(csb-cs)-rho*rs-(cs/L)*(Yac*rs*L+delta);
f7 = @(L,cco2b,cco2,rs,delta,Dco2) (Dco2/(Ll*L))*(cco2b-cco2)+4*rho*rs-(cco2/L)*(Yac*rs*L+delta) ;
f8 = @(L,chb,ch,rs,delta,Dh) (Dh/(Ll*L))*(chb-ch)+12*rho*rs-(ch/L)*(Yac*rs*L+delta) ;
f9 = @(rs,L,delta) -Am*(Yac*rs*L+delta) ;
f10 = @(csb,cs,vl,Ds) (1/vl)*(-(Am*Ds/Ll)*(csb-cs)) ;
f11 = @(cco2b,cco2,vl,Dco2) (1/vl)*(-(Am*Dco2/Ll)*(cco2b-cco2)) ;
f12 = @(i) -2.24*10^(-7)*(i)+0.0067;
%f12 = @(i) (kla*co2equi+7.7167*10^(-8)*exp(-(i)*(kla+qo2)))/(kla+qo2);
f13 = @(co2,ch)  Eocathode-((R*T*1000)/(4*F))*log(1/(co2*(ch^4))) ;
f14 = @(ch,cco2,cs) Eoanode-((R*T*1000)/(12*F))*log(((cco2^3)*(ch)^(11))/cs) ;
f15 = @(csb,Ds) (12*F*Ds*1000*csb)/(Ll) ;
f16 = @(I) ((dm/km)+(dcell/kaq))*(I/1000) ;
f17 = @(I,il) ((R*T)/(12*F))*1000*log(il/(il-I)) ;
f18 = @(I,cs) (b/2.303)*asinh(I/(2*ioref*cs)) ;
f19 = @(Ecathode,Eanode,nohm,nconc,nact) (abs(Ecathode-Eanode)-nohm-nconc-nact)*0.081 ;
f20 = @(Eoutput)  Eoutput/100;


d_t=1;
t = 0:d_t:18000;

mu=zeros(1,length(t));
kd=zeros(1,length(t));
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
ovr=zeros(1,length(t));
idenstity=zeros(1,length(t));
OUR=zeros(1,length(t));
xbyxo=zeros(1,length(t));
T=zeros(1,length(t));
power=zeros(1,length(t));
mutnoexp=zeros(1,length(t));
eff=zeros(1,length(t));

mu(1)=1.166370409*10^(-3);
rs(1)=2.5*10^(-20);%INITIAL I NEED TO BE KNOWN(TAKEN ZERO HERE) to vary
phia(1)=0.4286802857;%to vary
L(1)=0.0258*10^(-1);%to vary
cs(1)=7.5;
delta(1)=-1.446759259*10^(-08);
cco2(1)=1.804718175*10^(-10);%INITIAL CO2 CONC NEEDED and here up needed to be reviewed
ch(1)=10^(-7);
vl(1)=135*10^(-6);
csb(1)=7.5;
cco2b(1)=1.804718175*10^(-25);
co2b(1)=7.23*10^(-3);
co2(1)=7.23*10^(-3);
Ecathode(1)=845.8563841;
Eanode(1)=845.8563841;
il(1)=8519.75;
nohm(1)=0;
nconc(1)=0;
nact(1)=0;
Eoutput(1)=0;
I(1)=0;
T(1)=308;
Ds(1) = (0.2544*10^(-4))/86400;
Dh(1) = (2.3328*10^(-4))/86400;
Dco2(1) = (0.9960*10^(-4))/86400;



% for j=1:(length(mL))
%      L(1)= mL(j);
%       T = mT(j)
    for i=1:(length(t)-1)
        
        DsEff(i) = (1-phia(i))*Ds(i);
        DhEff(i) = (1-phia(i))*Dh(i);
        Dco2Eff(i) = (1-phia(i))*Dco2(i);
        phia(i+1) = phia(i) + 0.1*d_t*f3(rs(i),phia(i),delta(i),L(i));
        L(i+1) = L(i) + 0.008*d_t*f4(rs(i),L(i),delta(i));
        delta(i+1) = delta(i) + 0.8*d_t*f5(L(i));
        cs(i+1) = cs(i) + 0.00031*d_t*f6(L(i),csb(i),cs(i),rs(i),delta(i),DsEff(i));
        cco2(i+1) = cco2(i) + d_t*f7(L(i),cco2b(i),cco2(i),rs(i),delta(i),Dco2Eff(i));
        ch(i+1) = ch(i) + 0.000000003*d_t*f8(L(i),chb,ch(i),rs(i),delta(i),DhEff(i));
        vl(i+1) = vl(i) + d_t*f9(rs(i),L(i),delta(i));
        csb(i+1) = csb(i) + 0.0003*d_t*f10(csb(i),cs(i),vl(i),DsEff(i));
        cco2b(i+1) = cco2b(i) + d_t*f11(cco2b(i),cco2(i),vl(i),Dco2Eff(i));
        co2(i+1) = f12(i+d_t);
        Ecathode(i+1) = f13(co2(i+1),10^(-7));
        Eanode(i+1) = f14(ch(i+1),cco2(i+1),cs(i+1));
        il(i+1) = f15(csb(i+1),DsEff(i));
        nohm(i+1) = f16(I(i));
        nconc(i+1) = f17(I(i),il);
        nact(i+1) = f18(I(i),cs(i+1));
        mu(i+1) = f1(cs(i+1));
        rs(i+1) = f2(mu(i+1),nact(i+1),phia(i+1));
        Eoutput(i+1) = f19(Ecathode(i+1),Eanode(i+1),nohm(i+1),nconc(i+1),nact(i+1));
        I(i+1) = f20(Eoutput(i+1));
        Ds(i+1)= DsEff(i);
        Dh(i+1)= DhEff(i); 
        Dco2(i+1) = Dco2Eff(i); 
        if i<=16000
        kd(i)=mu(1)*exp(-(16000/t(i)));
        elseif i>16000
            kd(i)=mu(1);
        end
        ph = -log10(ch);

    end
    mcs(j,:)=cs; 
    mph(j,:)=ph;
    mL_cal(j,:)=L;
    ovr=(nact+nconc+nohm);
    movr(j,:)=ovr;
    mporosity(j,:)=1-phia;
end 
ovr=(nact+nconc+nohm);
OUR=co2*qo2;

idensity=I/54*10^(-1);
mut=(mu(1)-kd).*t;
mutnoexp(:,:)=exp(mut);
power=Eoutput.*I;
eff(:,:)=power/76.7;
%T=308-717.56r^2;
%plot(idenstity(100:8000),ovr(100:8000))

%plot(cs(100:8000),ovr(100:8000))

%plot(I,ovr)
%plot(cs,ovr)
%plot(cs(1:5000),I(1:5000))
% plot(cs(1:12000),eff(1:12000))
%plot(cs,mu)
%plot(cs(1:16000),mut(1:16000))
% plot(cs(1:16000),mut(1:16000))
%plot(t(1:12000),power(1:12000))
%plot(cs(1:5000),OUR(1:5000))
%plot(t(1:10000),co2(1:10000))


figure (1)
% plot(t(1:14000),Eoutput(1:14000))
plot(ph(1:14000),cs(1:14000),'k')
xlabel('pH in anode biofilm','FontWeight','bold');
ylabel('Concentration of lactate (kg/m^{3})','FontWeight','bold')
xlim([3,7]);
ylim([2.5,7.55]);
saveas(gcf,'pH.tiff')

figure (2)
plot(ovr(1:14000),cs(1:14000),'k')
xlabel('Overpotential (mV)','FontWeight','bold');
ylabel('Concentration of lactate (kg/m^{3})','FontWeight','bold')
xlim([0 140]);
ylim([2.5,7.55]);
saveas(gcf,'overpotential.tiff')

% figure (3)
% % plot(cs(1:14000),eff(1:14000))
% plot(cs(1:12000),eff(1:12000))
% figure (4)
% plot(cs(1:14000),I(1:14000))

figure (5)
plot(L(1:14000),cs(1:14000),'k')
xlabel('Biofilm Thickness (m)','FontWeight','bold');
ylabel('Concentration of lactate (kg/m^{3})','FontWeight','bold')
xlim([2.6*10^(-3),3.3*10^(-3)]);
saveas(gcf,'biofilm_thickness.tiff')
% ylim([2,7.5]);

figure (6)
plot(mut(1:16000),cs(1:16000),'k')
xlabel('ln(C/C_o)','FontWeight','bold');
ylabel('Concentration of lactate (kg/m^{3})','FontWeight','bold')
xlim([0 11.8]);
ylim([0,8]);
saveas(gcf,'growth_kinetics_bact.tiff')

figure(7)
plot(mut(1:14500),L(1:14500));
xlabel('ln(C/C_o)','FontWeight','bold');
ylabel('Biofilm Thickness (m)','FontWeight','bold');
saveas(gcf,'growth_bact_With_thickness.tiff')

figure (8)
plot(L(1:14000),ch(1:14000),'k')
xlabel('Biofilm Thickness (m)','FontWeight','bold');
ylabel('Concentration of hydrogen ions (kg/m^{3})','FontWeight','bold')
xlim([2.6*10^(-3),3.3*10^(-3)]);
saveas(gcf,'biofilm_thickness_hydrogen.tiff')

figure (9)
plot(L(1:14000),cco2(1:14000),'k')
xlabel('Biofilm Thickness (m)','FontWeight','bold');
ylabel('Concentration of carbon dioxide biofilm (kg/m^{3})','FontWeight','bold')
xlim([2.6*10^(-3),3.3*10^(-3)]);
saveas(gcf,'biofilm_thickness_CarbonDioxide.tiff')

figure (10)
plot(t(1:5401),cs(1:5401),'k')
xlabel('time(sec)','FontWeight','bold');
ylabel('Concentration of lactate (kg/m^{3})','FontWeight','bold')
xlim([0,5400]);
% ylim([2.5,7.55]);
saveas(gcf,'timevsConc90min.tiff')


% 
% %%%%%%%generating the surface plot for anode potential vs current vs
% %%%%%%%substrate%%%%%
% % [X Y] = meshgrid(Eanode(1:14000),cs(1:14000))
% 
% % figure(1)
% % % surf(Eanode(1:14000),cs(1:14000),I(1:14000))
% % scatter3(Eanode(1:14000),cs(1:14000),I(1:14000))
% % xlabel('anode potential(mV)');
% % ylabel('Concentration of lactate (kg/m^{3})')

% figure (7)
% plot(Eanode(11:14000),I(11:14000),'k','LineWidth',0.5)
% xlabel('anode potential(mV)','FontWeight','bold')
% ylabel('Current(mA)','FontWeight','bold')
% saveas(gcf,'current_AnodePotential.tiff')


% ae=[0 0.00000000001 0.00000000002 0.00000000004 0.00000000006 0.00000000008 0.0000000001 0.0000000002 0.0000000004 0.0000000006 0.0000000008  0.000000001 0.000000002 0.000000004 0.000000006 0.000000008 0.00000001 0.00000002 0.00000004 0.00000006 0.00000008 0.0000001 0.0000002 0.0000004 0.0000006 0.0000008 0.000001 0.000002 0.000004 0.000006 0.000008 0.00001 0.00002 0.00004 0.00006 0.00008 0.0001 0.0002 0.0003 0.0004 0.0005 0.0006 0.0007 0.0008 0.0009 0.001];
% mx= [25.7232 25.8523 25.9580 25.7913 26.0469 26.1236 26.1910 26.4411 26.7420 26.9349 27.0772 27.1902 27.5516 27.9239 28.1452 28.3034 28.4266 28.8124 29.2019 29.4312 29.5946 29.7217 30.1184 30.5181 30.7532 30.9206 31.0508 31.4571 31.8662 32.1068 32.2781 32.4114 32.8271 33.2456 33.4917 33.6669 33.8032 34.2283  34.4783 34.6563 34.7946 34.9079 35.0038 35.0870 35.1605 35.2263];
% plot(ae, mx,'k','LineWidth',0.5)
% xlabel('Cathode Aeration Rate(kg/s)','FontWeight','bold')
% ylabel('Maximum Power(MicroWatts)','FontWeight','bold')
% saveas(gcf,'aerationrate_power.tiff')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for bulk analysis%%%%%%%%%%%%%%%%%%%%
% figure (1)
% % plot(t(1:14000),Eoutput(1:14000))
% plot(ph(1:14000),csb(1:14000),'k')
% xlabel('pH in anode biofilm','FontWeight','bold');
% ylabel('Concentration of lactate in bulk (kg/m^{3})','FontWeight','bold')
% xlim([3,7]);
% ylim([2.5,7.55]);
% saveas(gcf,'pH_bulk.tiff')
% figure (2)
% plot(ovr(1:14000),csb(1:14000),'k')
% xlabel('Overpotential (mV)','FontWeight','bold');
% ylabel('Concentration of lactate in bulk(kg/m^{3})','FontWeight','bold')
% xlim([0 140]);
% ylim([2.5,7.55]);
% saveas(gcf,'overpotential_bulk.tiff')
% % figure (3)
% % % plot(cs(1:14000),eff(1:14000))
% % plot(cs(1:12000),eff(1:12000))
% % figure (4)
% % plot(cs(1:14000),I(1:14000))
% figure (5)
% plot(L(1:14000),csb(1:14000),'k')
% xlabel('Biofilm Thickness (m)','FontWeight','bold');
% ylabel('Concentration of lactate bulk(kg/m^{3})','FontWeight','bold')
% xlim([2.6*10^(-3),3.3*10^(-3)]);
% saveas(gcf,'biofilm_thickness_bulk.tiff')
% % ylim([2,7.5]);
% figure (6)
% plot(mut(1:16000),csb(1:16000),'k')
% xlabel('ln(C/C_o)','FontWeight','bold');
% ylabel('Concentration of lactate bulk (kg/m^{3})','FontWeight','bold')
% xlim([0 11]);
% ylim([0,8]);
% saveas(gcf,'growth_kinetics_bact_bulk.tiff')
% 
% plot(t(1:14000),csb(1:14000))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%with multiple temperatures%%%%%%%%%%%%%%%%


% figure (1)
% % plot(t(1:14000),Eoutput(1:14000))
% plot(mph(1,1:14000),mcs(1,1:14000),'-k',mph(2,1:14000),mcs(2,1:14000),'--r',mph(3,1:14000),mcs(3,1:14000),'-.')
% xlabel('pH in anode biofilm','FontWeight','bold');
% ylabel('Concentration of lactate (kg/m^{3})','FontWeight','bold')
% xlim([3,7]);
% ylim([2.5,7.55]);
% legend('0.0025','0.004','0.001');
% legend('300K','303K','306K');
% saveas(gcf,'mpH_Temp.tiff')
% figure (2)
% plot(movr(1,1:14000),mcs(1,1:14000),'-k',movr(2,1:14000),mcs(2,1:14000),'--r',movr(3,1:14000),mcs(3,1:14000),'-.')
% xlabel('Overpotential (mV)','FontWeight','bold');
% ylabel('Concentration of lactate (kg/m^{3})','FontWeight','bold')
% xlim([0 140]);
% ylim([2.5,7.55]);
% % legend('0.0025','0.004','0.001');
%  legend('300K','303K','306K');
% saveas(gcf,'moverpotential_temp.tiff')
% % figure (3)
% % % plot(cs(1:14000),eff(1:14000))
% % plot(cs(1:12000),eff(1:12000))
% % figure (4)
% % plot(cs(1:14000),I(1:14000))
% figure (5)
% plot(mL_cal(1,1:14000),mcs(1,1:14000),'-k',mL_cal(2,1:14000),mcs(2,1:14000),'--r',mL_cal(3,1:14000),mcs(3,1:14000),'-.')
% xlabel('Biofilm Thickness (m)','FontWeight','bold');
% ylabel('Concentration of lactate (kg/m^{3})','FontWeight','bold')
% xlim([2.6*10^(-3),3.3*10^(-3)]);
% ylim([2.3,8]);
% legend('300K','303K','306K');
% saveas(gcf,'mbiofilm_thickness_temp.tiff')
% 
% figure (6)
% plot(mporosity(1,1:14000),mcs(1,1:14000),'-k',mporosity(2,1:14000),mcs(2,1:14000),'--r',mporosity(3,1:14000),mcs(3,1:14000),'-.')
% xlabel('Biofilm porosity','FontWeight','bold');
% ylabel('Concentration of lactate (kg/m^{3})','FontWeight','bold')
% xlim([0.0001,0.571319714300000]);
% % ylim([2.3,8]);
% legend('300K','303K','306K');
% saveas(gcf,'mbiofilm_porosity_temp.tiff')
% ylim([2,7.5]);

% figure (6)
% plot(mut(1:16000),cs(1:16000),'k')
% xlabel('ln(C/C_o)','FontWeight','bold');
% ylabel('Concentration of lactate (kg/m^{3})','FontWeight','bold')
% xlim([0 11.8]);
% ylim([0,8]);
% saveas(gcf,'growth_kinetics_bact.tiff')
% figure(7)
% plot(mut(1:14000),L(1:14000));
% xlabel('ln(C/C_o)','FontWeight','bold');
% ylabel('Biofilm Thickness (m)','FontWeight','bold');
% saveas(gcf,'growth_bact_With_thickness.tiff')
% figure (8)
% plot(L(1:14000),ch(1:14000),'k')
% xlabel('Biofilm Thickness (m)','FontWeight','bold');
% ylabel('Concentration of hydrogen ions (kg/m^{3})','FontWeight','bold')
% xlim([2.6*10^(-3),3.3*10^(-3)]);
% saveas(gcf,'biofilm_thickness_hydrogen.tiff')
% figure (9)
% plot(L(1:14000),cco2(1:14000),'k')
% xlabel('Biofilm Thickness (m)','FontWeight','bold');
% ylabel('Concentration of carbon dioxide biofilm (kg/m^{3})','FontWeight','bold')
% xlim([2.6*10^(-3),3.3*10^(-3)]);
% saveas(gcf,'biofilm_thickness_CarbonDioxide.tiff')
