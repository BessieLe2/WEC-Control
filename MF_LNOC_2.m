clear all
close all
clc
%Wave Radiation Model Parameters
A_r=[0 0 -17.9;1 0 -17.7;0 1 -4.41];
B_r=[38.6 379 89];
C_r=[0;0;1]';
D_r=0;
n_r=size(A_r,1);

%Prediction step
n_p=5;

step=1000;%simulation steps
%Wave information:Height h, velocity v, force f
%a significant wave height of 4 m, a peak period of 6 s, and a peakedness factor of 3.3.
Height=4;
T=6;
Gamma=3.3;
[h,v,f]=Wave_JONSWAP_with_F_H_V(Height,T,Gamma);
%Stage cost parameter
r=5*1e-3;
t_s=0.2;
L=f.*v*t_s;

