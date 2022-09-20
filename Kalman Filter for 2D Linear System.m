%% Initialization
clc;
clear all;
close all;
load('sensor_data_kalman_2D.mat')

%% Parameters
% Sampling Time
T=0.01;

% State Transistion Matrix
F = [1 T T^2/2 0 0 0; 
     0 1 T 0 0 0;
     0 0 1 0 0 0;
     0 0 0 1 T T^2/2;
     0 0 0 0 1 T;
     0 0 0 0 0 1;];

% Output Transistion Matrix
H = [1 0 0 0 0 0;
     0 0 0 1 0 0;];

% Process Covariance
sig = 0.1;
Q = sig^2*[ T^4/4, T^3/3, T^2/2 0 0 0;
           T^3/3, 2*T^3, T^2 0 0 0;
           T^2/2, T^2, 100*T 0 0 0;                            
           0 0 0 T^4/4, T^3/3, T^2/2;
           0 0 0 T^3/3, 2*T^3, T^2;
           0 0 0 T^2/2, T^2, 100*T;];

% Measurement Covariance
r = 0.1;
R = [r 0;
     0 r;];

Ttime = 629;

% Initial States
x_curr = zeros(6,1); 
x_curr(1,1) = -0.1658;
x_curr(4,1) = 0.9285;
x_curr(3,1) = 0.107413865;
x_curr(6,1) = 0.913175337;
x_curr(2,1) = 0.5;
x_curr(5,1) = 0.36;

% Covariance Matrix
P = eye(6);

%% Estimation using Position data

for i = 1:Ttime
  % State Estimation
  x_prev = x_curr;
  x1 = F*x_prev;
  % Covariance Estimation
  P_prev = Q +F*P*F';
  % Kalman Gain
  K = P_prev*H'*inv(H*P_prev*H'+R);
  % State Update
  measured_position = meas_pos(i,:);
  x_meas = measured_position';
  x_curr = x1+ K*(x_meas - H*x1);
  % Covariance Update
  P = P_prev - K*H*P_prev;
  x_est(1:6,i) = x_curr(1:6);
end 

x_1 = x_est(1,:)';
xv_1 = x_est(2,:)';
xa_1 = x_est(3,:)';
y_1 = x_est(4,:)';
yv_1 = x_est(5,:)';
ya_1 = x_est(6,:)';

% Plots
figure(1)
plot(1:629, x_true(:,1), '-b')
hold on
plot(1:629, x_1, '.-r');
title('X Position (true & estimated) vs Time');
xlabel('Time Step');
ylabel('X Position');
legend('True', 'Estimate');

figure(2)
plot(1:629, x_true(:,2), '-b');
hold on
plot(1:629, y_1,'.-r');
title('Y Position (true & estimated) vs Time');
xlabel('Time Step');
ylabel('Y Position');
legend({'True', 'Estimate'},'Location','best');

figure(3)
plot(1:629, xv_1, '.-r');
title('X Velocity vs Time');
xlabel('Time Step');
ylabel('X Velocity');


figure(4)
plot(1:629, yv_1, '.-r');
title('Y Velocity vs Time');
xlabel('Time Step');
ylabel('Y Velocity');

figure(5)
plot(1:629, xa_1, '.-r');
title('X Acceleration vs Time');
xlabel('Time Step');
ylabel('X Acceleration');

figure(6)
plot(1:629, ya_1, '.-r');
title('Y Acceleration vs Time');
xlabel('Time Step');
ylabel('Y Acceleration');

figure(13)
plot( x_1, y_1, '.-r');
title('X Position vs Y position');
xlabel('X Position');
ylabel('Y Position');

%% Estimation Using Acceleration Sensor
C = [0 0 1 0 0 0;
     0 0 0 0 0 1;];
x_curr = zeros(6,1); 
x_curr(1,1) = -0.1658;
x_curr(4,1) = 0.9285;
x_curr(3,1) = 0.107413865;
x_curr(6,1) = 0.913175337;
x_curr(2,1) = 0.5;
x_curr(5,1) = 0.36;

p= eye(6);

for i=1:Ttime
  % State Estimation
  x_prev = x_curr;
  x1 = F*x_prev;
  % Covariance Estimate
  P_prev = Q +F*P*F';
  % Kalman Gain
  K = P_prev*C'*inv(C*P_prev*C'+R);
  measured_acceleration = meas_acc(i,:);
  % State Update
  x_acc = measured_acceleration';
  x_curr = x1+K*(x_acc - C*x1);
  % Covariance Update
  P = P_prev - K*C*P_prev;
  x_est(1:6,i) = x_curr(1:6);
end 

x_2 = x_est(1,:)';
xv_2 = x_est(2,:)';
xa_2 = x_est(3,:)';
y_2 = x_est(4,:)';
yv_2 = x_est(5,:)';
ya_2 = x_est(6,:)';

% Plots
figure(7)
plot(1:629, x_2, '.-r');
hold on
plot(1:629, x_true(:,1), '-b')
title('X Position vs Time');
xlabel('Time Step');
ylabel('X Position');
legend({'True', 'Estimate'},'Location','best');

figure(8)
plot(1:629, y_2,'.-r');
hold on
plot(1:629, x_true(:,2), '-b')
title('Y Position vs Time');
xlabel('Time Step');
ylabel('Y Position');
legend('True', 'Estimate');

figure(9)
plot(1:629, xv_2, '.-r');
title('X Velocity vs Time');
xlabel('Time Step');
ylabel('X Velocity');

figure(10)
plot(1:629, yv_2, '.-r');
title('Y Velocity vs Time');
xlabel('Time Step');
ylabel('Y Velocity');

figure(11)
plot(1:629, xa_2, '.-r');
title('X Acceleration vs Time');
xlabel('Time Step');
ylabel('X Acceleration');

figure(12)
plot(1:629, ya_2, '.-r');
title('Y Acceleration vs Time');
xlabel('Time Step');
ylabel('Y Acceleration');

figure(13)
plot(x_2, y_2, '.-r');
title('X Postion vs Y Position');
xlabel('X Postion');
ylabel('Y Position');

figure(14)
plot(x_1, y_1, '.-r');
hold on
plot(x_2, y_2, '.-b');
title('XY Position from Pos Sensor and Acc Sensor');
xlabel('X Postion');
ylabel('Y Position');
legend({'Position Sensor', 'Acceleration Sensor'},'Location','best');
