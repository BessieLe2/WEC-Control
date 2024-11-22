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
Wave=Wave_JONSWAP(3,5,3.3);

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
max_iter=1000;
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

