% This script is a discrete-time implementation of a Kalman filter for a 
% second order linear oscillator with control and disturbance forces.  
% The notation follows that of OPTIMAL CONTROL AND ESTIMATION by 
% R. F. Stengel.
clc
clear
close all

%% Continuous time system parameters

H = @(X) [1    sin(X(1)).*tan(X(2))    cos(X(1)).*tan(X(2));...
          0               cos(X(1))              -sin(X(1));...
          0    sin(X(1)).*sec(X(2))    cos(X(1)).*sec(X(2))];

f = @(X,U) H(X)*U;
l = @(X,w_tilda) H(X)*w_tilda;

h = @(X) [-sin(X(2));sin(X(1)).*cos(X(2));X(3)];


%% Jacobian

A = @(X,U) [U(2)*cos(X(1))*tan(X(2))-U(3)*sin(X(1))*tan(X(2))                  U(2)*sin(X(1))*(sec(X(2))^2)+U(3)*cos(X(1))*(sec(X(2))^2) 0; ...
                               -U(2)*sin(X(1))-U(3)*cos(X(1))                                                                          0 0; ...
            U(2)*cos(X(1))*sec(X(2))-U(3)*sin(X(1))*sec(X(2))  U(2)*sin(X(1))*(sec(X(2))*tan(X(2)))+U(3)*cos(X(1))*(sec(X(2))^tan(X(2))) 0];

C = @(X) [                  0            -cos(X(2)) 0; ...
          cos(X(1))*cos(X(2))  -sin(X(1))*sin(X(2)) 0; ...
                            0                     0 1];
L = @(X) H(X);


%% Noise Characteristics

sigma_W0 = [10 10 10]'.*.3*pi/180;
sigma_V = [0.4 0.4 0.3]'*0.4;

W0 = diag(sigma_W0.^2);
V = diag(sigma_V.^2);

OMEGA = pi;
Tfinal = 10*pi/OMEGA;
dt = 0.001;
t = 0:dt:Tfinal;
U = [cos(OMEGA.*t);sin(OMEGA.*t);t.*0];

w = diag(sigma_W0)*randn(3,length(t));

% sigma_v = V^0.5;
v = diag(sigma_V)*randn(3,length(t));


sigma_X = [5 5 5]'.*pi/180;
P0 = diag(sigma_X.^2);
X0 = [10 20 14]'.*pi/180;
Xh = [0 0 0]';
X  = X0;
Y  = h(X0)+v(1);
P = P0;

% sigma_w = W0^0.5;


for n = 2:length(t)
    % Actual State
    X(:,n) = X(:,n-1) + dt.* H(X(:,n-1))*(U(:,n-1)+w(:,n-1));
    Y(:,n) = h(X(:,n)) + v(:,n);
    % step 1
    Xm = Xh(:,n-1) + dt.* (f(Xh(:,n-1),U(:,n-1)));
    % step 2
    Pm = P + dt.* (A(Xh(:,n-1),U(:,n-1))*P + P*A(Xh(:,n-1),U(:,n-1))' + L(Xh(:,n-1))*W0*L(Xh(:,n-1))');
    % step 3
    Gk = Pm*C(Xm)'/(C(Xm)*Pm*C(Xm)+V);
    % step 4
    Xh(:,n) = Xm + Gk * (Y(:,n)-h(Xm));
    % step 5
    P = (eye(3) - Gk*C(Xh(:,n)))*Pm*(eye(3) - Gk*C(Xh(:,n)))+Gk*V*Gk';
end
%% Plots
close all
c = 180/pi;
figure;
subplot(3,1,1);plot(t,X(1,:)*c,t,Xh(1,:)*c);legend({'$x_1$','$\hat{x}_1$'},'Interpreter','Latex')
ylabel('\phi (deg)')
title('States (With Process Noise)');
subplot(3,1,2);plot(t,X(2,:)*c,t,Xh(2,:)*c);legend({'$x_2$','$\hat{x}_2$'},'Interpreter','Latex')
ylabel('\theta (deg)')
subplot(3,1,3);plot(t,X(3,:)*c,t,Xh(3,:)*c);legend({'$x_3$','$\hat{x}_3$'},'Interpreter','Latex')
xlabel('Time (s)');
ylabel('\psi (deg)')

figure;
subplot(3,1,1);plot(t,U(1,:)*c+w(1,:)*c,t,U(1,:)*c);legend({'$p+\tilde{p}$','$p$'},'Interpreter','Latex')
ylabel('p (deg/s)')
title('Control (With Process Noise)');
subplot(3,1,2);plot(t,U(2,:)*c+w(2,:)*c,t,U(2,:)*c);legend({'$q+\tilde{q}$','$q$'},'Interpreter','Latex')
ylabel('q (deg/s)')
subplot(3,1,3);plot(t,U(3,:)*c+w(3,:)*c,t,U(3,:)*c);legend({'$r+\tilde{r}$','$r$'},'Interpreter','Latex')
xlabel('Time (s)');
ylabel('r (deg/s)')

figure;
subplot(3,1,1);plot(t,Y(1,:),t,Y(1,:)-v(1,:));legend({'$y_1$','$h_1(X)$'},'Interpreter','Latex')
ylabel('y_1 (g)')
title('Outputs and Measurement Noise (With Process Noise)');
subplot(3,1,2);plot(t,Y(2,:),t,Y(2,:)-v(2,:));legend({'$y_2$','$h_2(X)$'},'Interpreter','Latex')
ylabel('y_2 (g)')
subplot(3,1,3);plot(t,Y(3,:)*c,t,Y(3,:)*c-v(3,:)*c);legend({'$y_3$','$h_3(X)$'},'Interpreter','Latex')
xlabel('Time (s)');
ylabel('y_3 (deg)')

