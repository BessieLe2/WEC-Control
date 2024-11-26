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
%rank(ctrb(A_init,B_u))
%Augmented system parameters
n_x=size(A_init,1);
n_p=2;%prediction length
D=zeros(1,n_p);
D(1)=1;
I=ones(n_p-1);
T=[zeros(n_p-1,1),I;0,zeros(1,n_p-1)];
A=[A_init,B_w*D;zeros(n_p,n_x),T];
B=[B_u;zeros(n_p,1)];
%rank(ctrb(A, B)) 
%The wave excitation force is generated using a JONSWAP spectrum [45], with a significant wave high of 3 m,
% a peak period of 5 s, and a peakedness parameter of 3.3.
Height=3;
T_wave=5;
Gamma=3.3;
Wave=Wave_JONSWAP(Height,T_wave,Gamma);

%Stage cost parameters
%r=0.0011;
%t_s=0.1;
%L=uC_xX+0.5Ru^2
%R=2*t_s*r+2*C_z*B_u;
R = 0.0011;
C_X=[C_z*(A_init-eye(size(A_init))),C_z*B_w*D];


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
tolerance=1e-6;
max_iter=1000;
j=1;%Iterative order number
l=n_x+n_p+1;
%%%Initialize RLS parameters
n_theta=l*(l+1)/2;%%n_theta is the number of RLS parameter 
%Initialize input and output in RLS
Z=zeros(n_theta,num_steps);%Initialize input Z
Y=zeros(num_steps,1);%Initialize output Y
N = 120;                % Number of iterations (warming-up period)
lambda = 0.98;          % Forgetting factor
P = eye(n_theta);       % Initial covariance matrix (identity matrix)
Theta = zeros(n_theta, 1); % Initial parameter estimates

% Recursive Least Squares Implementation
Theta_history = zeros(N, n_theta); % Store parameter estimates
error_history = zeros(N, 1);       % Store errors
while norm(M-M0)>tolerance &&j<max_iter
    M0=M;
    F0=F;
    %Generate random u at the first training epoch
    if j==1
        u_random=200*randn(num_steps,1);
    end
%Generate output Y and input Z for a iteration
for k=1:num_steps-1
    X(:,k)=[x(:,k);Wave(k:k+n_p-1)'];
    x(:,k+1)=A_init*x(:,k)+B_u*u(:,k)+B_w*Wave(k);
    X_upper(:,k+1)=[x(:,k+1);Wave(k+1:k+n_p-1)';0];
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
   Y(k,:)=2*u_temp*C_X*X(:,k)+R*u_temp.^2;
end
 
%fprintf("Output: %s\n", mat2str(Y));

%RLS function
for k = 2:N
    % Create regression vector Z_k
    Z_k = Z(:,k);
    
    % Update P_k+1
    K_k = (P * Z_k) / (lambda + Z_k' * P * Z_k); % Kalman gain
    P = (P / lambda) - K_k * Z_k' * P / lambda;
    
    % Update Theta_k+1
   
    prediction_error = Y(k) - Z_k' * Theta;
    Theta = Theta + P * Z_k * prediction_error;
    if Theta
      OO=1;
    else
        disp(P)
    end

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

Muu=M(l,l)
MuX=M(l,1:l-1);
F=-inv(Muu)*MuX;
%Update iterative order number
j=j+1;
end


%Method 1 :use command idare to solve DARE
%[H,-F,L]=idare(A,B,0,R,C_X',eye(n_x+n_p));



%Method 2
tolerance = 1e-6;                      % Convergence threshold
max_iter = 1000;                       % Maximum number of iterations
% Initialization
H = zeros(n_x+n_p);                   % Initial guess for H
H_prev = H + 2;                       % Ensure the loop starts
iter = 0;
 
% Iterative computation
while norm(H - H_prev, 'fro') > tolerance && iter < max_iter
     H_prev = H;                        % Save previous H
     G = C_X + B' * H_prev * A;         % Gain term
     S = R + B' * H_prev * B;           % Modified cost
     H = A' * H_prev * A - G' * (S \ G); % Update H
     iter = iter + 1;                   % Increment iteration counter
end
F=-inv(R+B'*H*B)*(C_X+B'*H*A);
 % Check for convergence
if iter >= max_iter
     disp('Warning: Iterative solution did not converge.');
 else
     disp(['Converged in ', num2str(iter), ' iterations.','F=',num2str(F)]);
 end
% 
% % Output solution
% disp('Solution H:');
% disp(H);
%Method 3
%Define initial F and H
H=eye(n_x+n_p);
H=[-2.4677    0.0006    0.0001    0.0004; 0.0006   -0.2078   -0.0001   -0.0001;0.0001   -0.0001   -0.0000   -0.0000;0.0004   -0.0001   -0.0000   -0.0000]*1e3;
H_prev=H+2;
F=[10 -1 3 5];
%F=[82.2145  -66.0540    0.0153    0.0549];

F_prev=F+2;
%F=[82.2149  -66.0544    0.0153    0.0549];
tolerance=1e-4;
max_iter=1000;
iter=1;
while norm(H-H_prev,'fro')>tolerance ||norm(F-F_prev,'fro')>tolerance &&iter<max_iter
    H_prev=H;
    F_prev=F;
    H=(A+B*F_prev)'*H_prev*(A+B*F_prev)+F_prev'*R*F_prev+2*F_prev'*C_X;
    F=-inv(R+B'*H*B)*(C_X+B'*H*A);
    iter=iter+1;
end













