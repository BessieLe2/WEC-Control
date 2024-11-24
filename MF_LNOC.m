clear all
close all
clc
 
% MF LNOC
%%Simulation Case 1-With Full State Information
%Reduced-order model parameters
A_init=[0.7726 0.1834;-2.1783 0.7614];
B_u=[0.0588;0.5635]*10^(-3);
B_w=B_u;
C_z=[1,0];

rank(ctrb(A_init,B_u))
%Augmented system parameters
n_x=size(A_init,1);
n_p=2;%prediction length
D=zeros(1,n_p);
D(1)=1;
I=ones(n_p-1);
T=[zeros(n_p-1,1),I;0,zeros(1,n_p-1)];
A=[A_init,B_w*D;zeros(n_p,n_x),T];
B=[B_u;zeros(n_p,1)];
rank(ctrb(A, B)) 
%The wave excitation force is generated using a JONSWAP spectrum [45], with a significant wave high of 3 m,
% a peak period of 5 s, and a peakedness parameter of 3.3.
Height=3;
T_wave=5;
Gamma=3.3;
Wave=Wave_JONSWAP(Height,T_wave,Gamma);

%Stage cost parameters
r=0.0011;
t_s=0.1;
%L=uC_xX+0.5Ru^2
R=2*t_s*r+2*C_z*B_u;
C_X=[C_z*(A_init-ones(size(A_init))),C_z*B_w*D];

% [H,F,L]=idare(A,B,0,R,C_X',eye(n_x+n_p))% couldn't be solved


% 
% 
% tolerance = 1e-6;                      % Convergence threshold
% max_iter = 10000000;                       % Maximum number of iterations
% 
% % Initialization
% H = zeros(n_x+n_p);                   % Initial guess for H
% H_prev = H + 2;                       % Ensure the loop starts
% iter = 0;
% 
% % Iterative computation
% while norm(H - H_prev, 'fro') > tolerance && iter < max_iter
%     H_prev = H;                        % Save previous H
%     G = C_X + B' * H_prev * A;         % Gain term
%     S = R + B' * H_prev * B;           % Modified cost
%     H = A' * H_prev * A - G' * (S \ G); % Update H
%     iter = iter + 1;                   % Increment iteration counter
% end
% 
% % Check for convergence
% if iter >= max_iter
%     disp('Warning: Iterative solution did not converge.');
% else
%     disp(['Converged in ', num2str(iter), ' iterations.']);
% end
% 
% % Output solution
% disp('Solution H:');
% disp(H);
%Define initial F and H
H=zeros(n_x+n_p);
%F=[0,-50,0 0];
F=[81.2804 -65.2976 0.0148 0.0537];
F_prev=[0 0 0 0];
tolerance=1e-4;
max_iter=10;
H_prev=H+2;
iter=0;
while norm(H-H_prev,'fro')>tolerance&&norm(F-F_prev,'fro')>tolerance&&iter<max_iter
    H_prev=H;
    H=(A+B*F)'*H_prev*(A+B*F)+F'*R*F+2*F'*C_X;
    F_prev=F;
    F=-inv(R+B'*H*B)*(C_X+B'*H*A);
    iter=iter+1;
end
disp("error of H:"+norm(H-H_prev,'fro'));
disp("error of F:"+norm(F-F_prev,'fro'));

%Generate stage cost by wave information
% Constants
rho = 1025;    % Density of seawater (kg/m^3)
g = 9.81;      % Gravitational acceleration (m/s^2)

% Calculate wave energy density (J/m^2)
L = (1/8) * rho * g * Wave.^2;
Y=2*L;

num_steps=1000;%Simulations steps
n_u=1;%Input number
%Initialize xk, Xk and X_upper k+1
X=zeros(n_x+n_p,num_steps);
X_upper=zeros(n_x+n_p,num_steps+1);
x=zeros(n_x,num_steps);
u=zeros(n_u,num_steps);
%Initialize M0 and F0
M=zeros(n_x+n_p+1);
F=zeros(1,n_x+n_p);
M0=M+1;
F0=F+1;

j=1;%Iterative order number
l=n_x+n_p+1;
%%%Initialize RLS parameters
n_theta=l*(l+1)/2;%%n_theta is the number of RLS parameter 
%Initialize input and output in RLS
Z=zeros(num_steps,n_theta);%Initialize input Z
Y=zeros(num_steps,1);%Initialize output Y
N = 120;                % Number of iterations (warming-up period)
lambda = 0.98;          % Forgetting factor
P = zeros(n_theta);       % Initial covariance matrix (identity matrix)
Theta = zeros(n_theta, 1); % Initial parameter estimates

% Recursive Least Squares Implementation
Theta_history = zeros(N, n_theta); % Store parameter estimates
error_history = zeros(N, 1);       % Store errors
while norm(F-F0)>tolerance | norm(M-M0)>tolerance &&j<max_iter


for k=1:num_steps-1
    X(:,k)=[x(:,k);Wave(k:k+n_p-1)'];
    x(:,k+1)=A_init*x(:,k)+B_u*u(:,k)+B_w*Wave(k);
    X_upper(:,k+1)=[x(:,k+1);Wave(k+1:k+n_p-1)';0];
    %Generate Z
    if j==1
        u_temp=200*randn(1,num_steps+1);
        z=[X(:,k);u_temp(k)];
        z_upper=[X_upper(:,k+1);u_temp(k+1)];
    else
       z=[X(:,k);F0*X(:,k)];
       z_upper=[X_upper(:,k+1);F0*X_upper(:,k+1)];
    end   
   
        %z=[X;FX],l=nx+np+1
        i=1;
        for p=1:l
            for q=p:l
                Z(k,i)=z(p)*z(q)-z_upper(p)*z_upper(q);
            end
        end
    %Generate Y
    u_temp=F0*X(:,k);
    Y(k)=2*u_temp*C_X*X(:,k)+R*u_temp^2;
end
for k = 2:N
    % Create regression vector Z_k
    Z_k = Z(k,:)';
    
    % Update P_k+1
    K_k = (P * Z_k) / (lambda + Z_k' * P * Z_k); % Kalman gain
    P = (P / lambda) - K_k * Z_k' * P / lambda;
    
    % Update Theta_k+1
    prediction_error = Y(k) - Z_k' * Theta;
    Theta = Theta + P * Z_k * prediction_error;

    % Store results
    Theta_history(k, :) = Theta';
    error_history(k) = prediction_error;
end

%Update M and F
j=j+1;
disp(norm(F-F0))
end






