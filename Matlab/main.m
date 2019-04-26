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

W = 10; % Process noise covariance, for wheel vertical speed, in (N/kg)^2
V = 1e-8; % Measurement noise covariance,
          % for distance between wheel and vehicle, in (m)^2

G = [0 0 0 1]';

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

%% Classical control

[b,a] = ss2tf(A,B,eye(4),zeros(1,4)');

U2X(1,:) = tf(b(1,:),a);
U2X(2,:) = tf(b(2,:),a);
U2X(3,:) = tf(b(3,:),a);
U2X(4,:) = tf(b(4,:),a);

[mag1,phase1,wout1] = bode(U2X(1,:));
[mag2,phase2,wout2] = bode(U2X(2,:));
[mag3,phase3,wout3] = bode(U2X(3,:));
[mag4,phase4,wout4] = bode(U2X(4,:));

figure; bode(U2X(1,:)); title('TF - u to x1')
figure; bode(U2X(2,:)); title('TF - u to x2')
figure; bode(U2X(3,:)); title('TF - u to x3')
figure; bode(U2X(4,:)); title('TF - u to x4')

figure; plot(roots(b(1,:)),'o'); hold on
plot(roots(a),'*')
title('TF - u to x1')
legend('zeros','poles');

figure; plot(roots(b(2,:)),'o'); hold on
plot(roots(a),'*')
title('TF - u to x2')
legend('zeros','poles');

figure; plot(roots(b(3,:)),'o'); hold on
plot(roots(a),'*')
title('TF - u to x3')
legend('zeros','poles');

figure; plot(roots(b(4,:)),'o'); hold on
plot(roots(a),'*')
title('TF - u to x4')
legend('zeros','poles');


