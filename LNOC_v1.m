clear;
close all;
clc;
%radiation force
A_r=[0,0,-17.9;1,0,-17.7;0,1,-4.41];
B_r=[36.5,394,75.1]';
C_r=[0,0,1];
n_r=size(A_r,1);
%wave excitation force
A_e=[0 0 0 0 -409;1 0 0 0 -459;0 1 0 0 -226; 0 0 1 0 -64.0;0 0 0 1 -9.96];
B_e=[1549866 -116380 24748 -644 19.3]';
C_e=[0 0 0 0 1];
n_e=size(A_e,1);

k_s=3866;

m_s=242;
m_a=83.5;
m=m_a+m_s;
%Build a model
%A_c=[0 1 0 0;-k_s/m 0 C_r/m -C_e/m;0 B_r A_r 0;0 0 0 A_e];
%B_wc=[0 0 0 B_e];
%B_uc=[0 1/m 0 0];
%C_z=[0 1 0];%n_r,n_e

%Define continuous system
A_c=[0 1 zeros(1,n_r) zeros(1,n_e);
    -k_s/m 0 C_r/m -C_e/m;
    zeros(n_r,1) B_r A_r zeros(n_r,n_e);
   zeros(n_e,n_r+2) A_e];
B_wc=[0; 0 ;zeros(n_r,1); B_e];
B_uc=[0; 1/m; zeros(n_r,1); zeros(n_e,1)];
C_z=zeros(1,2+n_e+n_r);
C_z(2)=1;
% Define the sampling period
T_s = 0.1; %  T = 0.1 seconds

%[A,B_u]=c2d(A_c,B_uc,T_s)
%%Zero-order holder
% Discretize the system
A = expm(A_c * T_s); % Calculate Ad using matrix exponential
% Calculate integral for Buc and Bwc discretization
integral_Ac = inv(A_c) * (expm(A_c * T_s) - eye(size(A_c))); % Integral of exp(Ac*t) dt from 0 to T
B_u = integral_Ac * B_uc; % Discretized input matrix for control input u
B_w = integral_Ac * B_wc; % Discretized input matrix for disturbance input w
%[A,B_u]=c2d(A_c,B_uc,T_s);
rank(ctrb(A,B_u))
% Simulation parameters
num_steps = 500; % Number of simulation steps, Time=50s
n=size(A,1);% Number of state
n_u=1;%Number of control input
x = zeros(n, num_steps); % Initialize state vector
x_corrected=zeros(n,num_steps);
u=zeros(n_u,num_steps); % Control input 
w = 0.5*zeros(num_steps,1); % Disturbance input (constant or can be a sequence)
%%
%Stage cost 
Q=diag([6,9.8,zeros(1,n_r+n_e)]);
r=0.08;
%Observer weights
R_v=0.01;
R_w=1;
%Prediction time horizon
n_p=5;%10,20,30
t_p=n_p*T_s;

%V=A^TVA+Q-(B_u^TVA+C_z)^T(r+B_u^TVB_u)^{-1}(B_u^TVA+C_z)
[V,K1,L1,info]=idare(A,B_u,Q,r,C_z',eye(n));

K_x=-inv(r+B_u'*V*B_u)*(C_z+B_u'*V*A);
phi=(A+B_u*K_x)';

psi=[V*B_w];
for i=1:n_p-1
    psi=[psi, phi^(i)*V*B_w];
end
K_d=-inv(r+B_u'*V*B_u)*B_u'*psi;
%% term related s is ignored
%s=zeros(1,num_steps);
%for i=num_steps:n_p+1
%    s(:,i-1)=phi*s(:,i)+phi*V*w(i-1)
%end
K_s=-inv(r+B_u'*V*B_u)*B_u';
%s_=s(:,n_p);
w_p=w(1:n_p,:);
%%Kalman Filter
[P,K2,L2]=idare(A',C_z',B_w*R_w*B_w',R_v,zeros(n,1),eye(n));
L=-P*C_z'*inv(C_z*P*C_z'+R_v);
%%
x_corrected(:,1)=x(:,1);
y=0.5;%% observed state
% Run the simulation
for k = 1:num_steps-1
    % ignore K_s*(phi^(n_p-1))*s_
    %u(:,k)=K_x*x(:,k)+K_s*(phi^(n_p-1))*s_+K_d*w(:,0:n_p);
    %u(:,k)=K_x*x(:,k)+K_d*w_p;
    u(:,k)=K_x*x_corrected(:,k)+K_d*w_p;
    %State update
    x(:, k+1) = A * x(:, k) + B_u * u(:,k) + B_w * w(k,:);
    %Kalman Filter
    x_corrected(:,k+1)=x(:,k+1)+L*(y-C_z*x(:,k+1));
end
% Output calculation
z=C_z*x;
t_p=(0:n_p-1)*T_s;
EIG=zeros(n_p,1);
for i=0:n_p-1
    eigenvalues =eig(phi^i);
    max_eigenvalue = max(eigenvalues)
    EIG(i+1)=abs(max_eigenvalue);
end
K_s

figure
plot(t_p,EIG)


