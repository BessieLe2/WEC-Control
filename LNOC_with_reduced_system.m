clear all
close all
clc
A=[0.7726 0.1834;-2.1783 0.7614];
Bu=[0.0588;0.5635]*1e-3;
Bw=Bu;
Cz=[0 1];
%Stage cost 
Q=diag([6,9.8]);
r=0.08;
%Observer weights
Rv=0.01;
Rw=1;
%Prediction time horizon
T_s=0.1;
n_p=50;%10,20,30
t_p=n_p*T_s;
% Simulation parameters
num_steps = 500; % Number of simulation steps, Time=50s
n=size(A,1);% Number of state
n_u=1;%Number of control input
x = zeros(n, num_steps); % Initialize state vector
x_corrected=zeros(n,num_steps);
u=zeros(n_u,num_steps); % Control input 
%The wave excitation force is generated using a JONSWAP spectrum [45], with a significant wave high of 3 m,
% a peak period of 5 s, and a peakedness parameter of 3.3.
Height=2;
T_wave=8;
Gamma=3.3;
w=Wave_JONSWAP(Height,T_wave,Gamma);
%V=A^TVA+Q-(B_u^TVA+C_z)^T(r+B_u^TVB_u)^{-1}(B_u^TVA+C_z)
[V,K1,L1,info]=idare(A,Bu,Q,r,Cz',eye(n));
Kx=-K1;
phi=(A+Bu*Kx)';
psi=[V*Bw];
for i=1:n_p-1
    psi=[psi, phi^(i)*V*Bw];
end
Kd=-inv(r+Bu'*V*Bu)*Bu'*psi;
%% term related s is ignored
%s=zeros(1,num_steps);
%for i=num_steps:n_p+1
%    s(:,i-1)=phi*s(:,i)+phi*V*w(i-1)
%end
Ks=-inv(r+Bu'*V*Bu)*Bu';
%s_=s(:,n_p);
%%Kalman Filter
[P,K2,L2]=idare(A',Cz',Bw*Rw*Bw',Rv,zeros(n,1),eye(n));
L=-P*Cz'*inv(Cz*P*Cz'+Rv);
x_corrected(:,1)=x(:,1);
%%Assume observed state is velocity
% zz=A*sin(2*pi*f*t);
% v=A*2*pi*f*cos(2*pi*f*t);
t=0:0.1:100;
f=0.2;
y=2.5*2*pi*f*cos(2*pi*f*t);%% observed state
% Run the simulation
for k = 1:num_steps-1
    % ignore K_s*(phi^(n_p-1))*s_
    %u(:,k)=K_x*x(:,k)+K_s*(phi^(n_p-1))*s_+K_d*w(:,0:n_p);
    %u(:,k)=K_x*x(:,k)+K_d*w_p;
    u(:,k)=Kx*x_corrected(:,k)+Kd*w(k:k+n_p-1)';
    %State update
    x(:, k+1) = A * x(:, k) + Bu * u(:,k) + Bw * w(k);
    %Kalman Filter
    x_corrected(:,k+1)=x(:,k+1)+L*(y(k+1)-Cz*x(:,k+1));
end
% Output calculation
z=Cz*x;
t_p=(0:n_p-1)*T_s;
EIG=zeros(n_p,1);
for i=0:n_p-1
    eigenvalues =eig(phi^i);
    max_eigenvalue = max(eigenvalues)
    EIG(i+1)=abs(max_eigenvalue);
end

figure
plot(t_p,EIG)
for i=1:n_p-1
    Coffiecient =Ks*(phi^i);
    eig(Coffiecient)
end
