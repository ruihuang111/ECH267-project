clear; 
clc;
close all

Nx = 3;   %dimension of state variable
Ny = 1;   %dimension of output

Qn = 2.098;       %battery capacity
Ts = 1;           %sample period for discretization

%Fitted equivalent circuit parameter
Ro = 0.0781;    %internal resistance     (Omega)
R1 = 0.0490;    %polarization resistance (Omega)
C1 = 1445.5;    %polarization capacity   (F)

tao1 = R1*C1;     %time constant

%building state function: x(k+1)=Ax(k)+Bu, x=[Vp;SOC;Ro] 
A = [1-Ts/tao1 0 0; 0 1 0; 0 0 1];
B = [Ts/C1;-Ts/Qn/3600;0];

%Coefficient of Voc(SOC)
p = [18.0598  -54.9081   62.0223  -31.2780    6.7504   -0.1133    3.6656]; 
%Coefficient of dVoc(SOC)/dSOC
dp = polyder(p);

%Import NASA dataset
load('random_current_NASA.mat'); 
load('random_voltage_NASA.mat');
Ubt = random_voltage_NASA;
Ibt = random_current_NASA;

%Deriving reference SOC through ampere-hour method
SOC_Ah(1) = 1.0;
for k = 2:length(Ibt)
    SOC_Ah(k) = SOC_Ah(k-1)-Ts*Ibt(k)/Qn/3600; 
end

%Initializing the nonlinear observer
N = length(Ibt);
Xe = zeros(Nx,N); %store the estimated state variable of the ovserver
Ye = zeros(Ny,N); %store the estimated output of the ovserver

Xe(:,1) = [0; 1*0.5; Ro*0.5]; %initial guess of x_hat
Xk = Xe(:,1);                 %x_hat

%Observer gain
k1 = 0.126;
k2 = 0.144; 
k3 = 0.0063;
K = [k1*A(1,1) 0 0; 0 k2*A(2,2) 0; 0 0 k3*A(3,3)]; 

%Nonlinear observer
for k=2:N
    Dh = [-1 polyval(dp,Xk(2,1)) -Ibt(k)];               %h'
    Yk = polyval(p,Xk(2,1)) - Ibt(k)*Xk(3,1) - Xk(1,1);  %y_hat = Voc(SOC_hat) - I*Ro_hat - Vp_hat 
    Xk = A*Xk + B*Ibt(k) + K*Dh'*(Ubt(k) - Yk);          %x_hat(k+1) = A*x_hat(k) + B*u + K*h'^T*(y-y_hat)
    
    Xe(:,k) = Xk;     %store the estimated state variable x_hat(k)
    Ye(:,k) = Yk;     %store the estimated output y_hat(k)
    
end

%Extract desired variables
SOC_est=Xe(2,:);
R0_est=Xe(3,:);
U_est=Ye;

%Calculate mean error
SOC_error_mean=mean(abs(SOC_Ah-SOC_est));  
U_error_mean=mean(abs(Ubt'-U_est));        

%Plotting
plot(1:N,SOC_Ah,1:N,SOC_est)
figure
plot(1:N,Ro*ones(N,1),1:N,R0_est)
