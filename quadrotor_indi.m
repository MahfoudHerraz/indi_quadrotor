% ********************* quadrotor control with INDI ********************* %
close all; clearvars;

%% model and simulation parameters

s = tf('s');

% drone parameters
m = 4; % mass of drone
g = 9.8; % gravity constant
% drone inertia in body frame
Ix = 0.033; 
Iy = 0.033;
Iz = 0.066;
J = diag([Ix Iy Iz]);
l = 0.5; % length of drone arm
a =10^-6; % rotors thrust coefficient T = a*omega^2
b =2*10^-7; % rotors counter-torque coefficient t = (+-)b*omega^2

% simulation parameters
tf = 10; % final time
ts = 0.05; % sampling time
Nf = tf/ts+1; % simulation steps number
t = linspace(0,tf,Nf);
sigma = 0.1;

% control loop gains
Kw = diag([10,10,10]);
Kphi = diag([4,4,4]);

%% Initialization

% desired attitude
step_time = 1;
step_value = [pi/4;pi/4;pi/4];
phi_d = [zeros(3,step_time/ts) repmat(step_value,1,Nf-step_time/ts)];

% initial state and control
phi0 = [0 0 0]';
phi_old = phi0;
w_old = [0 0 0]';
w_dot_old = [0 0 0]';
u_old = [0 0 0]';
q_old = quaternion(phi_old','euler','XYZ','frame');

% variables history
phi_list = zeros(3,Nf);
phi_list(:,1) = phi_old;
w_list = zeros(3,Nf);
w_list(:,1) = w_old;
wd_list = zeros(3,Nf-1);
u_list = zeros(3,Nf);
u_list(:,1) = u_old;
omega_list = zeros(4,Nf);
omega_list(:,1) = [0;0;0;0];

%% Simulation loop

tic
for k = 1:Nf-1
    disp(k);
    % outer loop
    w_d = Kphi*(phi_d(:,k)-phi_old);
    wd_list(:,k) = w_d;
    % inner loop
    w_d_dot = Kw*(w_d-w_old);
    % indi control
    u = u_old + J*(w_d_dot-w_dot_old);
    u_old = u;
    u_list(:,k+1) = u;
    % motor rotation speeds
    Ga = [0 l*a 0 -l*a; l*a 0 -l*a 0; b -b b -b];
    omega_list(:,k+1) = sqrt(pinv(Ga)*u);

    % update angular acceleration
    delta = sigma*randn(3);
    w_dot_old = (eye(3)+delta)*w_d_dot;
    % update angular rates
    w_old = w_old + ts*w_dot_old;
    % update attitude euler angles with rotation matrix
    % rot_mat = [1 sin(phi_old(1))*tan(phi_old(2)) cos(phi_old(1)*tan(phi_old(2)));...
    %     0 cos(phi_old(1)) -sin(phi_old(1)); ...
    %     0 sin(phi_old(1))/cos(phi_old(2)) cos(phi_old(1))/cos(phi_old(2))];
    % phi_old = phi_old + ts*rot_mat*w_old;
    % update attitude euler angles with quaternions
    [q0,q1,q2,q3] = parts(q_old);
    Sq = -1/2*[q1 q2 q3; -q0 q3 -q2; -q3 -q0 q1; q2 -q1 -q0];
    q_dot = Sq*w_old;
    q_old = q_old + ts*quaternion(q_dot(1),q_dot(2),q_dot(3),q_dot(4));
    phi_old = euler(q_old,'XYZ','frame')';

    % attitude and angular rates history
    phi_list(:,k+1) = phi_old;
    w_list(:,k+1) = w_list(:,k)+ts*w_dot_old;
end
toc

%% Plot results

figure();
plot(t,phi_d(1,:),'--b','DisplayName','reference');
hold on
plot(t,phi_list(1,:),'-r','DisplayName','phi');
hold on
plot(t,phi_list(2,:),'-k','DisplayName','theta');
hold on
plot(t,phi_list(3,:),'-g','DisplayName','psi');
legend
title('attitude plot');
figure();
plot(t,u_list);
title('control torques');
figure();
plot(t,omega_list);
title('motors speeds');