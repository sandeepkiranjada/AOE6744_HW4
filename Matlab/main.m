%% HW4 Problem - 1

%% Model

M = 250; % Vehicle Mass in kgs
m = 10; % Wheel Mass in kgs

k1 = 150e3; % Passive Stiffness in N/m
b1 = 1e3;  % Passive Damping in N/(m/s)

k2 = 10e3; % Tire Stiffness in N/m
b2 = 0.1e3;  % Tire Damping in N/(m/s)

%% Noise characteristics

W = 10; % Process noise covariance, for wheel vertical speed, in (N/kg)^2
V = 1e-8; % Measurement noise covariance, for distance between wheel and vehicle, in (m)^2


%% State-Space Represetation

A = [0          0       1          0; ...
     0          0       0          1; ...
 -k1/M       k1/M   -b1/M       b1/M; ...
  k1/m -(k1+k2)/m    b1/m -(b1+b2)/m]; ...
  
B = [0  0  1/M  -1/m]';

C = [1 -1 0 0];

G = [0 0 0 1]';

