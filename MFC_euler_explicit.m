clear all 
close all 

T = 273:10:453;  
d_t=1;
t = 0:d_t:18000;
cs_temp = zeros(length(T),length(t));



parfor j = 1:length(T)
         
    cs_temp(j,:)=model(T(j));
    
end 

% plot(t,cs);

power=Eoutput.*I;

% PV=Eoutput.^2./vl*100;  power normalized by volume







