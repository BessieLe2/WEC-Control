clear all
close all
clc
%
% 
% MF LNOC
%%Simulation Case 1-With Full State Information
%Reduced-order model parameters
A_init=[0.7726 0.1834;-2.1783 0.7614];
B_u=[0.0588;0.5635]*10^(-3);
B_w=B_u;
C_z=[1,0];
%Augmented system parameters
n_x=size(A_init,1);
n_p=2;%prediction length
D=zeros(1,n_p);
D(1)=1;
I=ones(n_p-1);
T=[zeros(n_p-1,1),I;0,zeros(1,n_p-1)];
A=[A_init,B_w*D;zeros(n_p,n_x),T];
B=[B_u;zeros(n_p,1)];

%The wave excitation force is generated using a JONSWAP spectrum [45], with a significant wave high of 3 m,
% a peak period of 5 s, and a peakedness parameter of 3.3.
Wave=Wave_JONSWAP(3,5,3.3);

%Stage cost parameters
r=0.0011;
t_s=0.1;
%L=uC_xX+0.5Ru^2
R=2*t_s*r+2*C_z*B_u;
C_X=[C_z*(A_init-ones(size(A_init)));C_z*B_w*D];
[H,F,L]=idare(A,B,0,R,C_X')


