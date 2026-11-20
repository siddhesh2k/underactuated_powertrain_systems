% clear all;
% clc
% close all;
%% PMSM
J=0.035; %0.035 motor inertia
Ld=0.5e-4; % d-axis inductance
Lq=1.3e-4; % q-axis inductance
Rs=17e-3; % stator resis
P=8; %poles
b=11e-5; %resistance/friction on mechanical motor
psim=0.075; %0.075 magnetic flux of permanent magnet
TL=500;
zeta=4;
gamma=0.8;
% wnq=(1/(1-gamma))*Rs/Lq;
% wnd=(1/(1-gamma))*Rs/Ld;
% wns=(1/(1-gamma))*b/J;
wnd=2000;  %2.1286e+3;
wnq=wnd;
wns=11.7021;
% J=0.47e-4;
% Ld=7e-3; %d-axis inductance
% Lq=7e-3; %q-axis inductance
% Rs=2.98; %stator resistance
% P=2;
% TL=2; %load torque on mechanical
% b=11e-5; %resistance/friction on mechanical
% psim=0.125; %flux of permanent magnet

kcq=2*zeta*wnq*Lq-Rs;
tiq=(2*zeta*wnq*Lq-Rs)/(wnq^2*Lq);
kcd=2*zeta*wnd*Ld-Rs;
tid=(2*zeta*wnd*Ld-Rs)/(wnd^2*Lq);
kcs=(2*zeta*wns-b/J)/(1.5*P^2*psim/J);
tis=(2*zeta*wns-b/J)/wns^2;

kpvq=kcq;
kivq=kcq/tiq;
kpvd=kcd;
kivd=kcd/tid;
kp=kcs;
ki=kcs/tis;
Cs=50000; %tire sleep 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DriveLine
C_d=0.351;          % Drag coefficient   
A_f=2.323;          % Frontal area for drag calculation 
rho=1.28;           % Air density       [kg/m^3]

m_t=2242;           % Vehicle mass CVW  [kg]
R_w=0.35;           % Radius of Wheel   [m]

g_r=10.453;         % Gear Ratio        [-] 



%%%%% Tire Slip

CVW_Wf = 0.456;
CVW_Wt = 1-CVW_Wf;

L = 2.9;            % Wheelbase [m]
hg = 0.518;         % [m]
La = L * CVW_Wf;
Lb = L * CVW_Wt;     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Battery
A_bat = [9.7537055e-01, 4.504345e-03, -2.28211183e-01, 1.18833486e-01;
     -1.15646738e-05, 1.0000212e+00, -1.05093146e-04, 5.55260818e-05;
     9.82383686e-04, -1.8131455e-04, 1.00089836e+00, -4.73382618e-03;
     -4.72926431e-04, 8.7263179e-05, -4.3285588e-03, 1.00227789e+00];
B_bat = [3.64608411e-04;
     2.77601760e-05;
     -1.47637225e-05;
     7.10423501e-06];
C_bat = [1 0 0 0; 0 1 0 0];
D_bat = [0; 0];
Initial = [4.18464;-0.01;0.05710627;0.97185065];

%% Fuel Cell
%Parameters
w=logspace(-3,5);
T = 333; 
R = 8.314;
Pamb = 1;
va = 0.016;
vc = 0.0025;
ka = 0.004;
kc = 0.01;
fh2 =2;
pch2 = 0.99;
cfh2 = 6.85e-4;
fo2 = 20;
pco2 = 0.21;
cfo2 = 6.804e-4;
pch2o = 0.01;
Air_in = fo2*pco2*cfo2;
Fuel_in = fh2*pch2*cfh2;
h2o_in = fo2*pch2o*cfo2;
% Air_in = 0;
% Fuel_in = 0;
% h2o_in = 0;
c = 3.5;
f = 9.6485e+04;
N = 436; %no. of fuel cells
Rohmic = 0.14;
Ra = 10e-4;
I = 10; %200 is best
%268 surface area cm^2
%xdot = Ax+Bu+Fd; 


A = [-1/(Ra*c) 0 0 0; 0 -R*T*ka/va 0 0; 0 0 -R*T*kc/vc 0; 0 0 0 -R*T*kc/vc]; %state matrix
B = [0 0 0; R*T/va 0 0; 0 R*T/vc 0; 0 0 R*T/vc]; %input matrix comprising with respect to Air_in, Fuel_in, h20_in
F = [1/(268*c); -R*T/(2*va*f); -R*T/(4*vc*f); R*T/(2*vc*f)]; %disturbance load matrix with respect to I (current) 

%operating points
%u0 and d0 = 0
% x10 = 0;
% x20 = 0;
% x30 = 0;
% x40 = 0;

%u0 and d0 constants
x10 = I*Ra;
x20 = R*T*Fuel_in/ka - R*T*I/(ka*2*f);
x30 = R*T*Air_in/kc - R*T*I/(kc*4*f);
x40 = R*T*h2o_in/kc + R*T*I/(kc*4*f);


%y = Cx+Du+Ed;
C = [-N (R*T*N/(2*f*(x20+Pamb))) (R*T*N/(4*f*(x30+Pamb))) (-R*T*N/(2*f*(x40+Pamb)))];
D = [0 0 0];
E = [-N*Rohmic/268];

% Concatenate the disturbance matrices with the input matrices
B_aug = [B F];  % Augmented input matrix
D_aug = [D E];  % Augmented feedthrough matrix

sys = ss(A, B_aug, C, D_aug);
Gtf = tf(sys);
Gtf = minreal(Gtf);
Gp = Gtf(:,1:end-1);
Gd = Gtf(:,end);
%Gp(:,2:end) = 0;
%Gp = Gp(1,1);
%%
[Ag,Bg,Cg,Dg]=ssdata(Gp);
[Ag,Bg,Cg,Dg]=minreal(Ag,Bg,Cg,Dg);

[Adis,Bdis,Cdis,Ddis]=ssdata(Gd);
[Adis,Bdis,Cdis,Ddis]=minreal(Adis,Bdis,Cdis,Ddis);

%%
Mp=2;
Md=2;
%w0p=0.5;
%w0d=10;
wc=2; %2
w0p=wc;
w0d=10;  %5
Amp=1e-4;
Amd=1e-2;

%% define Wp (Performance filter)
kp=1;
Wptf=Sko_weight(Mp,w0p,Amp,1);
%Wptf=eye(2)*Wp1;
Wpss=ss(Wptf); 
[Ap,Bp,Cp,Dp]=ssdata(Wpss);
[Ap,Bp,Cp,Dp]=minreal(Ap,Bp,Cp,Dp);
%% define Wu (Actuator filter)
Wu1=tf(4/pi);
Wutf = eye(3)*Wu1;
%Wutf = Wu1;
Wuss=ss(Wutf);
[Au,Bu,Cu,Du]=ssdata(Wuss);
[Au,Bu,Cu,Du]=minreal(Au,Bu,Cu,Du);
%% define Wd (Robust uncertainty filter)
Wdtf=tf([1 w0d/Md],[Amd w0d]);
Wdss=ss(Wdtf);
[Ad,Bd,Cd,Dd]=ssdata(Wdss);
[Ad,Bd,Cd,Dd]=minreal(Ad,Bd,Cd,Dd);

%%
Paug = augw(Gp,Wptf,Wutf,Wdtf);
opts=hinfsynOptions('Display','on', 'Method','RIC');
[K,CL,Gamma,info]=hinfsyn(Paug,1,3,opts); %1,3
Ktf=tf(K);
[acp,bcp,ccp,dcp]=ssdata(K);
Ktf=minreal(Ktf);
S=inv(eye(1,1)+Gp*Ktf);
S=minreal(S,1e-3);
Ttf=eye(1)-S;
Ttf=minreal(Ttf,1e-3);
Y=Ktf*S;
Y=minreal(Y,1e-3);

%% Lookup table
  % Define motor parameters
T_max = 500; % Maximum (constant) torque in Nm
omega_base = 5800; % Base speed in rpm
omega_max = 21000; % Maximum speed in rpm

% Convert speeds to rad/s
omega_base_rad = omega_base * 2 * pi / 60;
omega_max_rad = omega_max * 2 * pi / 60;

% Define speed range
omega_2d = linspace(0, omega_max_rad, 1000);

% Calculate torque for each speed
Torque_2d = zeros(size(omega_2d));

% Constant torque region
for i = 1:length(omega_2d)
    if omega_2d(i) <= omega_base_rad
        Torque_2d(i) = T_max;
    else
        % Constant power region
        Torque_2d(i) = T_max * (omega_base_rad / omega_2d(i));
    end
end


% Save omega and T to workspace
omega_rad = omega_2d; % Rename for clarity
save('lookup_data.mat', 'omega_2d', 'Torque_2d');


% Define motor parameters
T_max = 500; % Maximum torque in Nm
omega_base = 5800; % Base speed in rpm
omega_max = 21000; % Maximum speed in rpm

% Convert speeds to rad/s
omega_base_rad = omega_base * 2 * pi / 60;
omega_max_rad = omega_max * 2 * pi / 60;

% Define speed range
omega = linspace(0, omega_max_rad, 1000);

% Initialize i_d reference array
i_d_ref = zeros(size(omega));

% Field weakening constant (example value, should be determined experimentally)
k_fw = -0.0001;

% Calculate i_d reference for each speed
for i = 1:length(omega)
    if omega(i) <= omega_base_rad
        i_d_ref(i) = 0; % No field weakening below base speed
    else
        i_d_ref(i) = k_fw * (omega(i) - omega_base_rad); % Linear field weakening
    end
end

% Save omega and i_d_ref to workspace
save('id_lookup_data.mat', 'omega', 'i_d_ref');

%load('KD_controller.mat')
%% Energy Management System (Hinf)
% System Parameters
% K_ems1 = 1;
% T_ems = 1;
% 
% % State-Space Model Definition
% A_ems = [1 K_ems1*T_ems; 0 1];
% B2_ems = [K_ems1*T_ems; 1];
% B1_ems = [K_ems1*T_ems; 0];
% C_ems = [0 -1];
% D2_ems = -1;
% D1_ems = 1;
% B_ems_aug = [B2_ems B1_ems];
% D_ems_aug = [D2_ems D1_ems];
% 
% % Define the sampling time
% Ts = 1; % or any other appropriate discrete sampling time
% 
% % Create State-Space Model
% sys_ems = ss(A_ems, B_ems_aug, C_ems, D_ems_aug,Ts);
% 
% % Convert to Transfer Function and Minimize
% Gtf_ems = d2c(tf(sys_ems));
% Gtf_ems = minreal(Gtf_ems);
% 
% % Separate Plant and Disturbance Models
% Gp_ems = Gtf_ems(:, 1:end-1);
% Gd_ems = Gtf_ems(:, end);
% 
% % Obtain State-Space Data
% [Ag_ems, Bg_ems, Cg_ems, Dg_ems] = ssdata(Gp_ems);
% [Ag_ems, Bg_ems, Cg_ems, Dg_ems] = minreal(Ag_ems, Bg_ems, Cg_ems, Dg_ems);
% [Adis_ems, Bdis_ems, Cdis_ems, Ddis_ems] = ssdata(Gd_ems);
% [Adis_ems, Bdis_ems, Cdis_ems, Ddis_ems] = minreal(Adis_ems, Bdis_ems, Cdis_ems, Ddis_ems);
% 
% %
% Mp_ems=2;
% Md_ems=2;
% %w0p=0.5;
% %w0d=10;
% wc_ems=200; %2
% w0p_ems=wc;
% w0d_ems=200;  %5
% Amp_ems=1e-4;
% Amd_ems=1e-4;
% % Define Weighting Functions
% Wptf_ems = Sko_weight(Mp_ems, w0p_ems, Amp_ems, 1);
% %Wptf_ems = c2d(ss(Wptf_ems), Ts); % Convert to discrete-time system with sampling time Ts
% 
% Wu1_ems = tf(4/pi);
% Wutf_ems = eye(1) * Wu1_ems;
% %Wutf_ems = c2d(ss(Wutf_ems), Ts); % Convert to discrete-time system with sampling time Ts
% 
% Wdtf_ems = tf([1 w0d_ems/Md_ems], [Amd_ems w0d_ems]);
% %Wdtf_ems = c2d(ss(Wdtf_ems), Ts); % Convert to discrete-time system with sampling time Ts
% 
% 
% % Augment the Plant with Weighting Functions
% Paug_ems = augw(Gp_ems, Wptf_ems, Wutf_ems, Wdtf_ems);
% 
% % H-infinity Synthesis
% opts = hinfsynOptions('Display', 'on', 'Method', 'RIC');
% [K_ems, CL_ems, Gamma_ems, info_ems] = hinfsyn(Paug_ems, 1, 1, opts);
% 
% % Convert to Transfer Function and Minimize
% Ktf_ems = tf(K_ems);
% Ktf_ems = minreal(Ktf_ems);
% [acp_ems,bcp_ems,ccp_ems,dcp_ems]=ssdata(K_ems);
% % Compute Sensitivity Functions
% S_ems = inv(eye(size(Gp_ems, 1)) + Gp_ems * Ktf_ems);
% S_ems = minreal(S_ems, 1e-3);
% Ttf_ems = eye(size(Gp_ems, 1)) - S_ems;
% Ttf_ems = minreal(Ttf_ems, 1e-3);
% Y_ems = Ktf_ems * S_ems;
% Y_ems = minreal(Y_ems, 1e-3);
% % sigma(Ttf_ems,w)
% % hold on
% % sigma(S_ems,w)
% % hold on
% % sigma(1/Wptf_ems,w)
% % hold on
% % sigma(1/Wdtf_ems,w)
% % legend('T','S','1/wp','1/wd')
% % sigma(Y_ems,w)
% 
% %save('simulink_pdem_output_stochastic_V10.mat', 'out')