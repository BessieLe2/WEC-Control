clear all
close all
clc
%Wave Radiation Model Parameters
A_r=[0 0 -17.9;1 0 -17.7;0 1 -4.41];
B_r=[38.6 379 89];
C_r=[0;0;1]';
D_r=0;
n_r=size(A_r,1);

%Wave information:Height h, velocity v, force f
%a significant wave height of 4 m, a peak period of 6 s, and a peakedness factor of 3.3.
Height=4;
T=6;
Gamma=3.3;
[H,Vel,Force]=Wave_JONSWAP_with_F_H_V(Height,T,Gamma);
%Stage cost parameter
r=5*1e-3;
t_s=0.2;

%%Heave elevation zk, and heave velocity vk, are assumed to be directly measurable.
A=2.5;
f=0.2;
t=0:0.1:100;
zz=A*sin(2*pi*f*t);
v=A*2*pi*f*cos(2*pi*f*t);

%%Simulation

n_p=5;%Prediction step
n_u=1;%Input number
n_x=2;%Reduced System's state number
num_steps=1000;%Simulations steps

%Initialize xk, Xk and X_upper k+1
%Initialize M0 and F0
M=zeros(n_x+n_p+1);
F=zeros(1,n_x+n_p);
M0=M+1;
F0=F+1;
tolerance=1e-6;
max_iter=1000;
j=1;%Iterative order number
l=n_x+n_p+1;

%%Initialize RLS parameters
n_theta=l*(l+1)/2;%%n_theta is the number of RLS parameter 
%Initialize input and output in RLS
Z=zeros(n_theta,num_steps);%Initialize input Z
Y=zeros(num_steps,1);%Initialize output Y
N = 1000;                % Number of iterations (warming-up period)
lambda = 0.98;          % Forgetting factor

%Initialize xk, Xk and X_upper k+1
X=zeros(n_x+n_p,num_steps);
X_upper=zeros(n_x+n_p,num_steps+1);
% Recursive Least Squares Implementation
Theta_history = zeros(N, n_theta); % Store parameter estimates
error_history = zeros(N, 1);       % Store errors
while (norm(M-M0,'fro')>tolerance||norm(F-F0,'fro')>tolerance)&&j<max_iter
    M0=M;
    F0=F;
    %Generate random u at the first training epoch
    if j==1
        u_random=200*randn(num_steps+1,1);
    end
%Generate output Y and input Z for a iteration
for k=1:num_steps
    X(:,k)=[zz(k);v(k);Force(k:k+n_p-1)'];
    X_upper(:,k+1)=[zz(k+1);v(k+1);Force(k+1:k+n_p-1)';0];
    %Generate Z
    if j==1
        u_temp=u_random(k);
        z=[X(:,k);u_random(k)];
        z_upper=[X_upper(:,k+1);u_random(k+1)];
    else
       u_temp=F0*X(:,k);
       z=[X(:,k);F0*X(:,k)];
       z_upper=[X_upper(:,k+1);F0*X_upper(:,k+1)];
    end 
    %z=[X;FX],l=nx+np+1
    i=1;
        for p=1:l
            for q=p:l
                Z(i,k)=z(p)*z(q)-z_upper(p)*z_upper(q);
                i=i+1;
            end
        end
   %Generat Y
   Y(k,:)=Force(k)*Vel(k)*t_s;
end
Z=Z/norm(Z);%Normalize output Z
 
P = eye(n_theta);       % Initial covariance matrix (identity matrix)
Theta = zeros(n_theta, 1); % Initial parameter estimates

%RLS function
for k = 1:N
    % Create regression vector Z_k
    Z_k = Z(:,k);
    
    % % Update P_k+1
    K_k = P * Z_k / (lambda + Z_k' * P * Z_k);  % Compute Kalman gain
    P = (P - K_k * Z_k' * P) / lambda;  % Update covariance matrix

    % Update Theta_k+1
    prediction_error = Y(k) - Z_k' * Theta;
    Theta =Theta + P * Z_k * prediction_error;
    
    % Store results
    Theta_history(k, :) = Theta';
    error_history(k) = prediction_error;
end

%Update M and F
order=0;
for p=1:l
    for q=p:l
        order=order+1;
        if p==q
            M(p,q)=Theta(order);
        else
            M(p,q)=0.5*Theta(order);
            M(q,p)=M(p,q);
        end
    end
end
Muu=M(l,l);
MuX=M(l,1:l-1);
F=-inv(Muu)*MuX;
%Update iterative order number
j=j+1;
end


