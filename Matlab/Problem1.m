%% HW4 Problem - 1

clear
clc
close all

%% Model

M = 250; % Vehicle Mass in kgs
m = 10; % Wheel Mass in kgs

k1 = 150e3; % Passive Stiffness in N/m
b1 = 1e3;  % Passive Damping in N/(m/s)

k2 = 10e3; % Tire Stiffness in N/m
b2 = 0.1e3;  % Tire Damping in N/(m/s)

%% Noise characteristics

W0 = 10; % Process noise covariance, for wheel vertical speed, in (N/kg)^2
V = 1e-8; % Measurement noise covariance,
          % for distance between wheel and vehicle, in (m)^2

Gamma = [0 0 0 1]';

W = Gamma*W0*Gamma';

%% State-Space Represetation

getA = @(Mf,mf,k1f,k2f,b1f,b2f) ... Returns 'A' Matrix given the parameters
         [     0             0         1             0 ; ...
               0             0         0             1 ; ...
         -k1f/Mf        k1f/Mf   -b1f/Mf        b1f/Mf ; ...
          k1f/mf -(k1f+k2f)/mf    b1f/mf -(b1f+b2f)/mf]; ...
  
getB = @(Mf,mf) [0  0  1/Mf  -1/mf]'; % Returns 'B' Matrix



A = getA(M,m,k1,k2,b1,b2);

B = getB(M,m);

C = [1 -1 0 0];

D = 0; % For completeness



%% Controllablity and Observablity

if (4==rank(ctrb(A,B)))
    disp('system is fully contollable');
else
    error('system is not fully contollable');
end
    
if (4==rank(obsv(A,C)))
    disp('system is fully observable');
else
    error('system is not fully observable');
end


%% sub-problem 1

close all


R3 = [0 0 1 0];

Q = diag([1e4 100 0 0]);
R = B(3,:)^2;
% R = 1;
rho = 1;
% rho = sqrt(R/B(3,:)^2);

A_bar = (eye(4) - B*R3./B(3,:))*A;
B_bar = B.*rho;

if (4==rank(ctrb(A_bar,B_bar)))
    disp('new system is fully contollable');
else
    error('new system is not fully contollable');
end
    
if (4==rank(obsv(A_bar,C)))
    disp('new system is fully observable');
else
    error('new system is not fully observable');
end


[K_bar,sk,ek] = lqr(A_bar,B_bar,Q,R);
K = K_bar.*rho-R3*A./B(3,:);
[G,sg,eg] = lqr(A',C',W,V);
G=G';

A_aug = [   A-B*K    B*K;...
         zeros(4)  A-G*C];
X0 = [0.1 0 0 0]';
X_hat_0 = [0 0 0 0]';
E0 = X0-X_hat_0;
X_aug_0 = [X0;E0];

tf = 10;
dt = 0.00001;
t=0:dt:tf;

sigma_v = V^0.5;
v = sigma_v*randn(length(t),1);

sigma_w = W0^0.5;
w = sigma_w*randn(length(t),1);

X_aug_d = X_aug_0;
X_aug_n = X_aug_0;

for n=2:length(t)
    Xdd = A_aug*X_aug_d(:,n-1);
    Xdn = Xdd+[Gamma.*w(n);Gamma.*w(n)-G.*v(n)];
    X_aug_d(:,n) = X_aug_d(:,n-1) + dt*Xdd;
    X_aug_n(:,n) = X_aug_n(:,n-1) + dt*Xdn;

end

Xd = X_aug_d(1:4,:);
X_hat_d = Xd - X_aug_d(5:8,:);
ud = -K*X_hat_d;
Xn = X_aug_n(1:4,:);
X_hat_n = Xn - X_aug_n(5:8,:);
un = -K*X_hat_n;
%%
close all
figure;
subplot(2,1,1);plot(t,Xd(1,:),t,X_hat_d(1,:));legend({'$x_1$','$\hat{x}_1$'},'Interpreter','Latex')
ylabel('State 1')
title('States 1 and 2, No disturbance');
subplot(2,1,2);plot(t,Xd(2,:),t,X_hat_d(2,:));legend({'$x_2$','$\hat{x}_2$'},'Interpreter','Latex')
xlabel('Time (s)');
ylabel('State 2')


figure;
subplot(2,1,1);plot(t,Xd(3,:),t,X_hat_d(3,:));legend({'$x_3$','$\hat{x}_3$'},'Interpreter','Latex')
ylabel('State 3')
title('States 3 and 4, No disturbance');
subplot(2,1,2);plot(t,Xd(4,:),t,X_hat_d(4,:));legend({'$x_4$','$\hat{x}_4$'},'Interpreter','Latex')
ylabel('State 4')
xlabel('Time (s)');


figure; plot(t,ud);legend({'$u$'},'Interpreter','Latex')
ylabel('Control')
xlabel('Time (s)');
title('Control input, No disturbance');

figure;
subplot(2,1,1);plot(t,Xn(1,:),t,X_hat_n(1,:));legend({'$x_1$','$\hat{x}_1$'},'Interpreter','Latex')
ylabel('State 1')
title('States 1 and 2, With disturbance');
subplot(2,1,2);plot(t,Xn(2,:),t,X_hat_n(2,:));legend({'$x_2$','$\hat{x}_2$'},'Interpreter','Latex')
xlabel('Time (s)');
ylabel('State 2')


figure;
subplot(2,1,1);plot(t,Xn(3,:),t,X_hat_n(3,:));legend({'$x_3$','$\hat{x}_3$'},'Interpreter','Latex')
ylabel('State 3')
title('States 3 and 4, With disturbance');
subplot(2,1,2);plot(t,Xn(4,:),t,X_hat_n(4,:));legend({'$x_4$','$\hat{x}_4$'},'Interpreter','Latex')
ylabel('State 4')
xlabel('Time (s)');


figure; 
subplot(2,1,1);plot(t,un);legend({'$u$'},'Interpreter','Latex')
ylabel('Control')
title('Control input and Noise, With disturbance');
subplot(2,1,2); plot(t,w,t,v);legend({'$w$','$v$'},'Interpreter','Latex')
ylabel('Noise, w an v')
xlabel('Time (s)');

%% Classical control
clear tf
pause
close all
[num_U2U,den_U2U] = ss2tf(A,B,K,0); % TF from U to U, LQR
UUg = tf(num_U2U,den_U2U); % Loop Gain fron U to U

[num_U2Y,den_U2Y] = ss2tf(A,B,C,0); % TF from U to Y = C(sI-A)^-1 B, LQG
UYg = tf(num_U2Y,den_U2Y); % Gain from U to Y

A_tilda = A-B*K-G*C;
[num_Y2Xh,den_Y2Xh] = ss2tf(A_tilda,G,K,0); % TF from Y to Xhat K(s) = K(sI-A-BK-GC)^-1 G, LQG
KofS = tf(num_Y2Xh,den_Y2Xh);

LQG = UYg*KofS;
LQR = UUg;

[G2,sg2,eg2] = lqr(A',C',W+1e1.*B*B',V);
G2=G2';

A_tilda2 = A-B*K-G2*C;

[num_Y2Xh,den_Y2Xh] = ss2tf(A_tilda2,G2,K,0); % TF from Y to Xhat
KofS2 = tf(num_Y2Xh,den_Y2Xh);
LQGLTR = UYg*KofS2;

[G2,sg2,eg2] = lqr(A',C',W+1e10.*B*B',V);
G2=G2';

A_tilda2 = A-B*K-G2*C;

[num_Y2Xh,den_Y2Xh] = ss2tf(A_tilda2,G2,K,0); % TF from Y to Xhat
KofS2 = tf(num_Y2Xh,den_Y2Xh);
LQGLTR2 = UYg*KofS2;

figure; bode(LQG); hold on;grid on
bode(LQR); bode(LQGLTR); legend show

[magq,phasq,wq] = bode(LQG);
[magr,phasr,wr] = bode(LQR);
[magqr,phasqr,wqr] = bode(LQGLTR);
[magqr2,phasqr2,wqr2] = bode(LQGLTR2);

wall = unique([wq;wr;wqr;wqr2]);

[magq,phasq,wq] = bode(LQG,wall);
[magr,phasr,wr] = bode(LQR,wall);
[magqr,phasqr,wqr] = bode(LQGLTR,wall);
[magqr2,phasqr2,wqr2] = bode(LQGLTR2,wall);

magq = reshape(magq,1,length(magq(1,1,:)))';
phasq = reshape(phasq,1,length(phasq(1,1,:)))';
magr = reshape(magr,1,length(magr(1,1,:)))';
phasr = reshape(phasr,1,length(phasr(1,1,:)))';
magqr = reshape(magqr,1,length(magqr(1,1,:)))';
phasqr = reshape(phasqr,1,length(phasqr(1,1,:)))';
magqr2 = reshape(magqr2,1,length(magqr2(1,1,:)))';
phasqr2 = reshape(phasqr2,1,length(phasqr2(1,1,:)))';

figure; 
subplot(2,1,1);semilogx(wq,mag2db(magq),'.-r'); hold on;grid on
semilogx(wr,mag2db(magr),'g');
semilogx(wqr,mag2db(magqr),'b');
semilogx(wqr,mag2db(magqr2),'k');
ylabel('Magnitude (dB)');
title('Bode plot')
axis tight
legend('LQG','LQR','LQG with LTR r=10','LQG with LTR r=1e10')
subplot(2,1,2);semilogx(wq,atan2(sind(phasq),cosd(phasq))*180/pi,'.-r'); hold on;grid on
semilogx(wr,atan2(sind(phasr),cosd(phasr))*180/pi,'g');
semilogx(wqr,atan2(sind(phasqr),cosd(phasqr))*180/pi,'b');
semilogx(wqr,atan2(sind(phasqr2),cosd(phasqr2))*180/pi,'k');
ylabel('Wrapped Phase (deg)');
xlabel('Frequency (rad/s)')
axis tight

figure; margin(LQG); grid on;legend('LQG')
figure; margin(LQR); grid on;legend('LQR')
figure; margin(LQGLTR); grid on;legend('LQG with LTR r=10')
figure; margin(LQGLTR2); grid on;legend('LQG with LTR r=1e10')
