% main model_mfc
clc
clear
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

constants();
y0 = initial_condition();
tend = 1000;
tspan = [0 tend];
M = M_matrix;
options = odeset('Mass',M,'RelTol',1);
[t,y]=ode23s(@model_mfc,tspan,y0,options);

