% This script is a discrete-time implementation of a Kalman filter for a 
% second order linear oscillator with control and disturbance forces.  
% The notation follows that of OPTIMAL CONTROL AND ESTIMATION by 
% R. F. Stengel.

clear
close all

% Continuous time system parameters
omegan = 1;
zeta = 0.1;
omegad = sqrt(1-zeta^2)*omegan;

A = [0, 1; -omegan^2, -2*zeta*omegan];
B = [0; omegan^2];
C = [1, 0];
L = [0; omegan^2];

% Choose the sample time to capture around Z points per cycle of damped
% oscillation.
Z = 10;
T = (1/Z)*((2*pi)/omegad);

% Discrete time system parameters
Phi = expm(A*T);
Gamma = (Phi - eye(2))*inv(A)*B;
Lambda = (Phi - eye(2))*inv(A)*L;

% Choose the maximum index so that the simulation captures N cycles of
% damped oscillation.
N = 10;
kmax = ceil((N/T)*(2*pi/omegad));

% Initialize the data arrays.
x = zeros(kmax+1,2);
xhat_minus = x;
xhat_plus = x;
y = zeros(kmax+1,1);
yhat = y;
P_minus = zeros(2,2,kmax+1);
P_plus = P_minus;
u = zeros(kmax+1,1);

% Gaussian random initial state
mean_x = 0;
sigma_x = 1;
x(1,:) = mean_x + sigma_x*randn(size(x(1,:)));
P_plus(:,:,1) = eye(2);

% Zero-mean random disturbance vector and (constant) covariance.
sigma_w = 10;
w = sigma_w*randn(kmax+1,1);
W = Lambda*(sigma_w^2)*Lambda';

% Zero-mean random noise vector and (constant) covariance.
sigma_v = 0.1;
v = sigma_v*randn(kmax+1,1);
V = sigma_v^2;

for k = 1:kmax
    % Actual state update
    x(k+1,:) = (Phi*x(k,:)' + Gamma*u(k) + Lambda*w(k))';
    % State estimate update based on dynamic model
    xhat_minus(k+1,:) = (Phi*xhat_plus(k,:)' + Gamma*u(k))';
    % State estimate covariance update based on dynamic model
    P_minus(:,:,k+1) = Phi*P_plus(:,:,k)*Phi' + W;
    % Measurement and measurement estimate
    y(k+1) = C*(x(k+1,:))' + v(k+1);
    yhat(k+1) = C*xhat_minus(k+1,:)';
    % Optimal observer gain
    G = P_minus(:,:,k)*C'*inv(C*P_minus(:,:,k)*C' + V);
    % State estimate update based on measurement
    xhat_plus(k+1,:) = (xhat_minus(k+1,:)' + G*(y(k+1)-yhat(k+1)))';
    % State estimate covariance update based on measurement
    P_plus(:,:,k+1) = inv(inv(P_minus(:,:,k+1)) + C'*inv(V)*C);
end

figure
title('State and State Estimates')
plot(T*[1:kmax+1]',x(1:kmax+1,1),'.',...
    T*[1:kmax+1]',xhat_minus(1:kmax+1,1),'-.',...
    T*[1:kmax+1]',xhat_plus(1:kmax+1,1),'-.','LineWidth',2)
legend('True State','Model-Updated Estimate','Measurement-Updated Estimate')
xlabel('Time (s)')

P11 = zeros(kmax+1,1);
P12 = zeros(kmax+1,1);
P22 = zeros(kmax+1,1);
for i = 1:kmax+1
P11(i) = P_plus(1,1,i);
P12(i) = P_plus(1,2,i);
P22(i) = P_plus(1,2,i);
end

% figure
% title('Error Covariance Matrix Elements')
% plot(T*[1:kmax+1]',P11,T*[1:kmax+1]',P12,T*[1:kmax+1]',P22)
% legend('P_{11}','P_{12}','P_{22}')
% xlabel('Time (s)')

% figure
% title('Output and Output Estimate')
% plot(T*[1:kmax+1]',y(1:kmax+1,1),'.',...
%     T*[1:kmax+1]',yhat(1:kmax+1,1),'--')
% legend('Output','Output Estimate')

