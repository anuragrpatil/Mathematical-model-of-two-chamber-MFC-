
 t=0:1:5*60*60;
 co2equi=7.26*10^(-3);
kla=414/86400;
qo2=2.64/86400;
Co2 = zeros(1,length(t));

for i=1:1:length(t)

Co2(i) = (kla*co2equi+7.7167*10^(-8)*exp(-t(i)*(kla+qo2)))/(kla+qo2);
%Co2(i) = -2.24*10^(-7)*(i)+0.0067;
end


