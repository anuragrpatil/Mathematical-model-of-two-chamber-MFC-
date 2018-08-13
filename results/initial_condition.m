function y0 = initial_condition()
y0 = zeros(20,1);
global bdet

y0(1)=0; %rs
y0(2)=0; %nu
y0(3)=0.0478*10^(-10); %phi_a
y0(4)=10^(-10); %phi_i
y0(5)=0.025; %L
y0(6)=-bdet*y0(5); %delta
y0(7)=5; %cs 
y0(8)=0.23*10^(-7); %cco2 
y0(9)=5; %ch 
y0(10)=135*10^(-6); %vl 
y0(11)=0; %csb 
y0(12)=0; %cco2b 
y0(13)=7.23*10^(-3); %co2 
y0(14)=0; %Eoutput 
y0(15)=0; %Ecathode 
y0(16)=0; %Eanode 
y0(17)=0; %nact
y0(18)=0; %nohm
y0(19)=0; %nconc
y0(20)=0; %il
y0(21)=0; %i
y0(22)=0; %chb