clc
clear
close all
%% Wave energy converter model:
ts = 0.2;
Ts = ts;
k_s=3866;
m_s=242;
m_a=83.5;
m=m_a+m_s;
Df = 20
%Define continuous system
A_c=[0 1 ;
    -k_s/m -Df/m ];
B_wc=[0; 0 ];
B_uc=[0; 1/m];
[A,B]=c2d(A_c,B_uc,Ts);
B_u=B;
B_w=B;
Cz=[1,0];

np=2;
D=zeros(1,np);
D(1)=1;
C_X=[C_z*(A-eye(size(A))),C_z*B_w*D];


%% read the wave data
load wave
% N=200;
% TS=0.5;
% Ts=ts;
% tend=500;
% trunc = 250/TS+1;
% time = wdata(1:trunc,1);
% w = wdata(2.5*trunc+1:trunc*3.5,2);
ti = (0:ts:Fex(end,1))'; 
wi=interp1(Fex(:,1),Fex(:,2),ti);
dim_wi = size(wi,1);
%%OR
%Wave information:Height h, velocity v, force f
%a significant wave height of 4 m, a peak period of 6 s, and a peakedness factor of 3.3.
Height=4;
T=6;
Gamma=3.3;
[H,Vel,Force]=Wave_JONSWAP_with_F_H_V(Height,T,Gamma);
%wi=Wave'

%% Stage cost parameter
R=1e-3;
Rloc = R + 2*Cz*B_u;

%% Augmented system
nx=2;
np=2;
nu=1;
l=nx+np+1;
ntheta=l*(l+1)/2;

%% Initialisating simulation
Nm = 7000;
xk=zeros(nx,1);
RecordX=zeros(nx,Nm);
RecordU=zeros(1,Nm);
RecordE=zeros(1,Nm);
RecordP=zeros(1,Nm);
ek = 0;
RecordError = zeros(1,Nm);
RecordH = zeros(ntheta,Nm);
thetaH = zeros(ntheta,1); % recursive LS
PthetaH = 0.99*eye(ntheta);
lambda= 0.98;
ZZ=zeros(ntheta,1);
M=zeros(l);
F= -[0 15 zeros(1,np)];
%F=[81.2804  -65.2976    0.0148    0.0537];
%F=[173.2 -136.9 -0.0147 -0.0552 0.0929 -0.0603 -0.0077];
%% staring runing onlines
for j=1:5
    for i=1:Nm
        Xk = [xk;wi(i:np+i-1)];
        uk = 50*rand();
        %Record current input and state vector
        xkm1 = xk;
        RecordX(:,i) = xk;
        RecordU(:,i) = uk;
        %Caculate state at next step
        xk=A*xk+B_u*uk+B_w*wi(i); % update states
        Xdkp1=[xk;wi(i+1:np+i-1);0];
      
        % record energy and power
        %     pk = uk*Cz*(xkm1-xk)-R*uk^2;
        pk = uk*Cz*(xkm1-xk)-0.5*R*uk^2;
        RecordP(:,i) = pk/ts;
        ek = ek +pk;
        RecordE(:,i) = ek;
        % Policy evaluation:
        z=[Xk;uk];
        z_upper=[Xdkp1;F*Xdkp1];
        k=1;
        for p=1:l
            for q=p:l
                ZZ(k)=z(p)*z(q)-z_upper(p)*z_upper(q);
                k=k+1;
            end
        end
    % (Zk-Zkm1)*thetaP = Lk  
    %NewY = Lk;
    %NewX = (Zk-Zkm1)';
    NewX=ZZ';
        Lk = -pk; % stage cost
    %ukp1k = F*Xdkp1;
    %Zkm1 = toZbar(Xk,uk);
    %Zk = toZbar(Xdkp1,F*Xdkp1);
    error=C_X*Xk+Cz*(xkm1-xk)
     Y=0.5*(2*uk*C_X*Xk+R*uk^2);
    % (Zk-Zkm1)*thetaP = Lk  
    NewY = Y;
    %NewX = (Zk-Zkm1)';

    %% policy evaluation using Recursive LS
    PthetaH = 1/lambda * PthetaH - 1/lambda * PthetaH*NewX'*inv(lambda+NewX*PthetaH*NewX')*NewX*PthetaH;
    thetaH = thetaH + PthetaH*NewX'*(NewY-NewX*thetaH);
    RecordH(:,i) = thetaH;
    RecordError(:,i) = NewY-NewX*thetaH;
    end
    figure()
    plot(RecordError)
    %  H = zeros(np+nu+nx);
    % nH = np+nu+nx;
    % Hvec = thetaH;        
    % pos = @(j)nH*(j)-(j)*(j-1)/2;
    % for j=1:nH  
    % 
    %     H(:,j) =  [zeros(j-1,1);Hvec(pos(j-1)+1:pos(j),1)];
    % end
    % H=(H'+H)/2;    
    % H=mat2cell(H,[nx+np,nu],[nx+np,nu]);
    % Huu=cell2mat(H(2,2));
    % Hux=cell2mat(H(2,1));
    % F=-inv(Huu)*Hux
    %KXRE = [KXRE;Kx]
    order=0;
    for p=1:l
        for q=p:l
            order=order+1;
            if p==q
                M(p,q)=thetaH(order);
            else
                M(p,q)=0.5*thetaH(order);
                M(q,p)=M(p,q);
            end
        end
    end
Muu=M(l,l);
MuX=M(l,1:l-1);
F=-inv(Muu)*MuX;
end