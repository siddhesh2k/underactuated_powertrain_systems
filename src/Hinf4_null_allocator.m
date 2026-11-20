clear;
close;
clc


%% Parameters
Q = 100; % Battery capacity in Ah
Em = 720; % Battery open circuit voltage
eta_fc = 0.5; % Fuel cell efficiency
Q_lhv = 120000; % Lower heating value of hydrogen in J/g

%% State Space
A_pwr1 = [0,0;0,0];
B_pwr1 = [1/(3600*Em*Q); 1/(eta_fc*Q_lhv)];
F_pwr1 = [-1/(3600*Em*Q); 0];
C_pwr1 = eye(2,2);
D_pwr1 = [0;0];
E_pwr1 = [0;0];



B_aug_pwr1 = [B_pwr1 F_pwr1];
D_aug_pwr1 = [D_pwr1 E_pwr1];

sys1 =ss(A_pwr1, B_aug_pwr1, C_pwr1, D_aug_pwr1);
Gtfpwr1 = tf(sys1);
Gdpwr1 = Gtfpwr1(:,end);
[Adpwr1,Bdpwr1,Cdpwr1,Ddpwr1]=ssdata(Gdpwr1);
num1 = 3.858e-9;
den1 = [1, 0]; 


num2 = 2.016 * [2.733, 1.115e06, 1.234e09, 3.211e11];
den2 = 720 * 1.45 * [1, 1.187e07, 1.948e10, 7.864e12, 3.515e13, 0]; 

G1 = tf(num1, den1);
G2 = tf(num2, den2);

Gppwr2 = [G1; G2];
[Appwr1,Bppwr1,Cppwr1,Dppwr1]=ssdata(Gppwr2);
Gpinv = inv(Gppwr2.'*Gppwr2)*Gppwr2.';
X1 = generateAllPass(Gpinv(1,1));
X2 = generateAllPass(Gpinv(1,2));
Gpinv(1,1) = makeproper(Gpinv(1,1),1e-2);
Gpinv(1,1) = minreal(Gpinv(1,1)*X1);
Gpinv(1,2) = minreal(Gpinv(1,2)*X2);
[Ainv, Binv, Cinv, Dinv] = ssdata(Gpinv);

%% weights 
% Wp1pwr_diag = makeweight(100000, [1e10, 1], 1/3);
% Wp1pwr = blkdiag(Wp1pwr_diag, Wp1pwr_diag);
% 
% Wutfpwr = makeweight(1/5,[1e10, 1/6], 1/10);%1e10
% 
% Wd1pwr_diag = makeweight(100000, [1e10, 1], 1/3);
% Wd1pwr = blkdiag(Wd1pwr_diag, Wd1pwr_diag);

Wp1pwr_diag = makeweight(1/5, [1e10, 1/6], 1/10);
Wp1pwr = blkdiag(Wp1pwr_diag, Wp1pwr_diag);

Wutfpwr = makeweight(1/5,[1e10, 1/6], 1/10);%1e10

Wd1pwr_diag = makeweight(1/10, [1e10, 1/6], 1/5);
Wd1pwr = blkdiag(Wd1pwr_diag, Wd1pwr_diag);

% Wp1pwr_diag = makeweight(1/1, [1e-2, 1/2], 1/3);
% Wp1pwr = blkdiag(Wp1pwr_diag, Wp1pwr_diag);
% 
% Wutfpwr = makeweight(1/5,[1e10, 1/6], 1/10);
% 
% Wd1pwr_diag = makeweight(1/1.5, [1, 2], 100);
% Wd1pwr = blkdiag(Wd1pwr_diag, Wd1pwr_diag);

% Wp1pwr_diag = makeweight(100, [1e-3, 1], 1/1.5);
% Wp1pwr = blkdiag(Wp1pwr_diag, Wp1pwr_diag);
% 
% Wutfpwr = makeweight(1/5,[1e10, 1/6], 1/10);
% 
% Wd1pwr_diag = makeweight(1/1.5, [1e-2, 2], 100);
% Wd1pwr = blkdiag(Wd1pwr_diag, Wd1pwr_diag);

% Wp1pwr_diag = makeweight(1/1, [1e1, 1/2], 1/3);
% Wp1pwr = blkdiag(Wp1pwr_diag, Wp1pwr_diag);
% 
% Wutfpwr = makeweight(1/5,[1e10, 1/6], 1/10);
% 
% Wd1pwr_diag = makeweight(1/1.5, [1, 2], 100);
% Wd1pwr = blkdiag(Wd1pwr_diag, Wd1pwr_diag);

% Wp1pwr_diag = makeweight(1/5, [1, 1/6], 1/10);
% Wp1pwr = blkdiag(Wp1pwr_diag, Wp1pwr_diag);
% 
% Wutfpwr = makeweight(1/5,[1, 1/6], 1/10);
% 
% Wd1pwr_diag = makeweight(1/10, [1, 1/6], 1/5);
% Wd1pwr = blkdiag(Wd1pwr_diag, Wd1pwr_diag);


%%
w=10^-5:0.01:10^3;
Paugpwr1 = minreal(augw(Gppwr2,Wp1pwr,Wutfpwr,Wd1pwr));
optspwr1=hinfsynOptions('Display','on', 'Method','RIC');
[Kpwr1,CLpwr1,Gammapwr1,infopwr1]=hinfsyn(Paugpwr1,2,1,optspwr1); 
Ktfpwr1=tf(Kpwr1);
[acppwr1,bcppwr1,ccppwr1,dcppwr1]=ssdata(Kpwr1);
Ktfpwr1=minreal(Ktfpwr1,1e-3);
Spwr1=inv(eye(2,2)+Gppwr2*Ktfpwr1);
%Spwr1=minreal(Spwr1,1e-3);
Tpwr1=eye(2)-Spwr1;
%Tpwr1=minreal(Tpwr1,1e-3);
Ypwr1=Ktfpwr1*Spwr1;
%Ypwr1=minreal(Ypwr1,1e-3);
figure(1)
bodemag(Tpwr1,w);
title('Ty')
% figure(1)
% sigma(Tpwr1)
% hold on
% sigma(Spwr1)
% hold on
% sigma(1/Wp1pwr)
% hold on 
% sigma(1/Wd1pwr)
% legend('T','S','1/wp','1/wd')
% hold on
% sigma(Ypwr1)
% % figure(3)
% % bodemag(Ktfpwr1)
% figure(2)
% sigma(Wp1pwr*Spwr1,w)
% hold on
% sigma(Wd1pwr*Tpwr1,w)
% hold on
% sigma(Wutfpwr*Ypwr1,w)
% legend('Wp*S','Wd*T','Wu*Y')
% figure(3)
% sigma(Spwr1*Gdpwr1)

%% 
% syms s
% %[UL, UR, SMF] = smithMcMillanForm_YH([3.858e-9; 1.389e-5], s);
% 
% %UL = sym2tf(UL);
% %UR = sym2tf(UR);
% UL = tf([1 0; -1.389e-5/3.858e-9 1]); % Did Smith mcmillan by hand



syms s 

Gp1=3.858e-9/s;
Gp2=720*2.016*(2.733*s^3 + 1.115e06*s^2 + 1.234e09*s + 3.211e11)/(1.45*(s^5 + 1.187e07*s^4+ 1.948e10*s^3 + 7.864e12*s^2 + 3.515e13*s));%1.389e-5/s
Gp=[Gp1;Gp2];

Gp=simplifyFraction(Gp);
[Num,Den] = numden(Gp);
P1 = simplifyFraction(Num);
del=lcm(Den);
P=simplifyFraction([P1(1,1)*(del/Den(1,1));P1(2,1)*(del/Den(2,1))]);

[UL, UR, SMF] = smithMcMillanForm(P, del);
UL = [0 1;1 0]*UL;
SMF = [0 1; 1 0]*SMF;
UL = sym2tf(UL);


%% Prefilter Design
%pre_filter design;
Ka2 = UL(1,:).'; %UL(2,:).';
W = [1 0;0 0];  %[1 -150;0 1]

% Gaa = Ka2*inv(inv(Ka2.'*Ka2)*Ka2.'*W*Ka2)*inv(Ka2.'*Ka2)*Ka2.'*W;
% Ga = Gaa; 
Tpwr1 = minreal(Tpwr1,1e-2);
Spwr1 = minreal(Spwr1,1e-2);

[Q1, ~] = schur(evalfr(Tpwr1, 0.02), 'real');  % Ensures real-valued decomposition

% Apply Schur transformation to Tpwr1
% Tpwr1_schur = Q1' * Tpwr1 * Q1;

Ka3 = Tpwr1*Ka2;
% Ka3(1,1) = makeproper(Ka3(1,1),1e-3); %remove when not using shur
% Ka3(2,1) = makeproper(Ka3(2,1),1e-3);  %remove when not using shur
Ga = minreal(Ka2*inv(Ka3.'*W*Ka3)*Ka3.'*W, 1e-1);
% H = tf( (2*pi*100)^2 , [1 sqrt(2)*2*pi*100 (2*pi*100)^2] )
H = tf( (2*pi*10)^2 , [1 1/sqrt(2)*2*pi*10 (2*pi*10)^2] );
% H = 10^(-3.81/20) * tf( (0.02)^2 , [1 1/sqrt(2)*0.02 (0.02)^2] );
allpass1 = generateAllPass(Ga(1,1));
allpass2 = generateAllPass(Ga(2,1));
% Ga(1,1) = makeproper(minreal(Ga(1,1)*H*allpass1,1e-1),1e-3);  %remove makeproper when not using shur
% Ga(2,1) = makeproper(minreal(Ga(2,1)*H*allpass2,1e-1),1e-3);  %remove makeproper when not using shur
Ga(1,1) = (minreal(Ga(1,1)*H*allpass1,1e-2));  %remove makeproper when not using shur
Ga(2,1) = (minreal(Ga(2,1)*H*allpass2,1e-2));  %remove makeproper when not using shur

[Aa,Ba,Ca,Da]=ssdata(Ga(1,1));
[Num1,Den1] = (tfdata(Ga(1,1)));
Num1 = Num1{1};
Den1 = Den1{1};
[Num2,Den2] = (tfdata(Ga(2,1)));
Num2 = Num2{1};
Den2 = Den2{1};
[Num3,Den3] = (tfdata(Tpwr1(2,1)*Ga(1,1)));
Num3 = Num3{1};
Den3 = Den3{1};
new = minreal(Tpwr1(1,2)*Tpwr1(2,1)*Ga(1,1)/Tpwr1(2,2),1e-1);
allpass3 = generateAllPass(new);
new = minreal(new*allpass3);
[Num4,Den4] = tfdata(new);
Num4 = Num4{1};
Den4 = Den4{1};
figure;
% bodemag(Ga(1,1));
% hold on
% bodemag(Tpwr1(1,1));
% hold on

bodemag(Tpwr1(1,1)*Ga(1,1))
hold on
grid on
bodemag(Tpwr1(2,1)*Ga(1,1))
xlim([1e-2 1e4])
legend('Ty(1,1)*Ga(1,1)','Ty(2,1)*Ga(1,1)')

%%

% Compute Singular Value Decomposition (SVD) of Tpwr1
% [U, Sigma, V] = svd(evalfr(Tpwr1, 0.01));  % Evaluate at s=0 for static structure
% 
% % Apply the transformation matrix V
% Tpwr1_svd = U' * Tpwr1 * V;
% 
% 
% % Compute the exact pre-compensator
% W = [1, 0;
%     -Tpwr1(2,1) / Tpwr1(1,1), 1];
% 
% W = [1, 0; Tpwr1(2,2)*Tpwr1(2,1)/(Tpwr1(1,2)*Tpwr1(2,1)-Tpwr1(1,1)*Tpwr1(2,2)) -Tpwr1(1,1)*Tpwr1(2,2)/(Tpwr1(1,2)*Tpwr1(2,1)-Tpwr1(1,1)*Tpwr1(2,2))];
% 
% 
% allpass4 = generateAllPass(W(2,1));
% allpass5 = generateAllPass(W(2,2));
% W(2,1) = minreal(W(2,1)*allpass4);
% W(2,2) = minreal(W(2,2)*allpass5);
% 
% % Apply W(s) to decouple the system
% Tpwr1_block = W * Tpwr1;
% 
% [Aw, Bw, Cw, Dw] = ssdata(W);
% figure;
% subplot(2,1,1);
% bodemag(Tpwr1(2,1), 'r'); hold on;
% bodemag(Tpwr1_schur(2,1), 'g--');
% bodemag(Tpwr1_svd(2,1), 'b--');
% bodemag(Tpwr1_block(2,1), 'm--');
% legend('Original Ty(2,1)', 'Schur-Ty(2,1)', 'SVD-Ty(2,1)', 'Block-Ty(2,1)');
% title('Bode Plot: Comparison of Ty(2,1) Decoupling');
% grid on;
% 
% subplot(2,1,2);
% bodemag(Tpwr1(2,2), 'r'); hold on;
% bodemag(Tpwr1_schur(2,2), 'g--');
% bodemag(Tpwr1_svd(2,2), 'b--');
% bodemag(Tpwr1_block(2,2), 'm--');
% legend('Original Ty(2,2)', 'Schur-Ty(2,2)', 'SVD-Ty(2,2)', 'Block-Ty(2,2)');
% title('Bode Plot: Comparison of Ty(2,2) Preservation');
% grid on;
% 
% figure;
% bodemag(Tpwr1_schur); hold on; bodemag(Tpwr1); hold on; bodemag(Ga)
% legend('Ty_schur', 'Ty','Ga')
% hold on
% bodemag(Tpwr1(1,1)*Ga(1,1))
% hold on
% bodemag(Tpwr1(1,2)*Ga(2,1))
% hold on
% bodemag(Tpwr1(2,1)*Ga(1,1))
% hold on
% bodemag(Tpwr1(2,2)*Ga(2,1))
% legend('Ty(1,1)*Ga(1,1)','Ty(1,2)*Ga(2,1)','Ty(2,1)*Ga(1,1)','Ty(2,2)*Ga(2,1)')
% figure(5)
% bodemag(Tpwr1(1,1)*Ga(1,2))
% hold on
% bodemag(Tpwr1(1,2)*Ga(2,2))
% hold on
% bodemag(Tpwr1(2,1)*Ga(1,2))
% hold on
% bodemag(Tpwr1(2,2)*Ga(2,2))
% legend('Ty(1,1)*Ga(1,2)','Ty(1,2)*Ga(2,2)','Ty(2,1)*Ga(1,2)','Ty(2,2)*Ga(2,2)')


%%

% alpha = - (Tpwr1(2,1) * Ga(1,1) ) / (Tpwr1(2,1) - (Tpwr1(1,1) * Tpwr1(2,2) / Tpwr1(1,2)));
% allpass4 = generateAllPass(alpha);
% 
% alpha = minreal(alpha*allpass4);
% alpha1 = minreal(alpha*allpass4,1e-2);
% [Num5,Den5] = tfdata(alpha1);
% Num5 = Num5{1};
% Den5 = Den5{1};
% isstable(alpha1)
% isproper(alpha1)
% beta = -Tpwr1(1,1)/Tpwr1(1,2)*alpha;
% allpass5 = generateAllPass(beta);
% beta = minreal(beta*allpass5);
% beta1 = minreal(beta*allpass5,1e-2);
% [Num6,Den6] = tfdata(beta1);
% Num6 = Num6{1};
% Den6 = Den6{1};
% isstable(beta1)
% isproper(beta1)
% Define frequency range for Bode analysis


% bodemag(alpha); hold on; bodemag(alpha1);
% bodemag(beta); hold on; bodemag(beta1);