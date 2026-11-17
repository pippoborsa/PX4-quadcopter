clear
close all
clc

%% system parameters 

g = 9.81;
m = 1.07;
L = 0.25*sin(pi/4);
k = 0.0183;
km = 0.00026;
h = 0.085;

I_roll_xx = 0.02;
I_pitch_yy = 0.02;
I_yaw_zz = 0.03;

I = diag([I_roll_xx I_pitch_yy I_yaw_zz]);

c = 0.05;                             % spherical hinge damping coefficient

C = diag([c c c]);

%% EoM building blocks

syms t
syms yaw pitch roll real
syms yaw_d pitch_d roll_d real
syms yaw_dd pitch_dd roll_dd real
syms u_yaw u_pitch u_roll real

q = [roll pitch yaw];
qd = [roll_d pitch_d yaw_d];
qdd = [roll_dd pitch_dd yaw_dd];

U = [u_roll, u_pitch, u_yaw];                                   % [PWM]

F = [k*(U(3)+U(2)-U(1));
    k*(U(3)-U(2)+U(1));
    k*(-U(3)+U(2)+U(1));
    k*(-U(3)-U(2)-U(1))
    ];

M = [km*(U(3)+U(2)-U(1));
    km*(U(3)-U(2)+U(1));
    km*(-U(3)+U(2)+U(1));
    km*(-U(3)-U(2)-U(1))
    ];

% we need to parametrize rotations: YAW Z -> PITCH Y' -> ROLL X''

iRo = RotAxisAngle([0, 0, 1], yaw)';
iiRi = RotAxisAngle([0, 1, 0], pitch)';
iiiRii = RotAxisAngle([1, 0, 0], roll)';

omega = yaw_d*iiiRii*iiRi*iRo*[0 0 1]'+...
        pitch_d*iiiRii*iiRi*[0 1 0]'+...
        roll_d*iiiRii*[1 0 0]';

bRo = iiiRii*iiRi*iRo; % passive rotation matrix from fixed frame to body frame!

jac_qd = jacobian(omega, qd);

T = 1/2*qd*jac_qd'*I*jac_qd*qd';

V = m*g*h*[0 0 1]*bRo*[0; 0; 1];

V_disp = matlabFunction(V, 'vars', [pitch, roll]);

% figure
% x = -2*pi:0.1:2*pi;
% y = -2*pi:0.1:2*pi;
% [X, Y] = meshgrid(x, y);
% Z = V_disp(X, Y);
% surf(X, Y, Z, 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.2);
% xlabel('pitch')
% ylabel('roll')
% zlabel('V')
% 
% figure;
% imagesc(x, y, Z);  % Display Z as an image with x and y as axes
% colorbar;          % Show color bar for reference
% axis xy;           % Keep axis orientation correct
% xlabel('pitch');
% ylabel('roll');
% colormap jet;

D = 1/2*qd*C*qd';

Q = [2*L*(F(2)+F(3)-F(1)-F(4));
     2*L*(F(1)+F(3)-F(2)-F(4));
     -M(3)-M(4)+M(1)+M(2)
     ];
%% EoM

%%%%=======================================================================
%%%%     Build Lagrange equations
%%%%     D(dT/dqd)/Dt - dT/dq + dV/dq - Q = 0 
%%%%=======================================================================

% Compute dT/dqd
dT_dqd = simplify(jacobian(T,qd))';

% Compute D(dT/dqd)/Dt
DdT_Dtdqd = simplify(jacobian(dT_dqd,t)) ...
            + simplify(jacobian(dT_dqd,q))  * qd' ...
            + simplify(jacobian(dT_dqd,qd)) * qdd';

% Compute dT/dq
dT_dq = simplify(jacobian(T,q))' ;  

% Compute dD/dqd
dD_dq = simplify(jacobian(D,qd))' ; 

% Compute dV/dq
dV_dq = simplify(jacobian(V,q))' ; 
       
% Set up equations
equations = simplify(DdT_Dtdqd - dT_dq + dD_dq + dV_dq - Q);

f = solve(equations, qdd);
state_space_equations = [qd'; f.roll_dd; f.pitch_dd; f.yaw_dd];
state_space_equations = simplify(state_space_equations);

A = jacobian(state_space_equations, [q, qd]); % jacobian w.r.t. state
B = jacobian(state_space_equations, U);       % jacobian w.r.t. control
A_lin = double(subs(A, [q, qd, U], zeros(1,9)));
B_lin = double(subs(B, [q, qd, U], zeros(1,9)));
C_lin = [diag([1 1 1]), diag([0 0 0])];
D_lin = zeros(3,3);
sys = ss(A_lin, B_lin, C_lin, D_lin);

% p = [-3 -4 -3 -4 -1 -2];
% K_pole_placement = place(A_lin, B_lin, p)
% 
Q = diag(ones(6,1));
R = diag(ones(3,1));
% K_LQR = lqr(A_lin, B_lin, Q, R)

% augmented system for full state integral action
% u = k*[roll pitch yaw roll_d pitch_d yaw_d roll_i pitch_i yaw_i];

A_aug = zeros(9,9);
A_aug(1:6, 1:6) = A_lin;
A_aug(7:9, 1:6) = C_lin;
B_aug = [B_lin; zeros(3,3)];
Q_aug = diag([100 100 100 1 1 1 10 10 10]);
R_aug = R.*0.001;

K_LQR_aug = lqr(A_aug, B_aug, Q_aug, R_aug)

model_to_simulink = matlabFunction(state_space_equations, 'file', 'quadcopter','vars', {q, qd, U});

A = A_lin;
B = B_lin;
C = eye(6);
D = D_lin;

Vd = diag([1000 1000 1000 10 10 10]);
Vn  = diag([0.1 0.1 0.1 0.1 0.1 0.1]);
[Kf, P, E] = lqe(A, Vd, C, Vd, Vn);
Ak = A-Kf*C;
Bk = [B, Kf];
Ck = eye(6);
Dk = 0*[B Kf];
kalman = ss(Ak, Bk, Ck, Dk);
kalmand = c2d(kalman, 0.01, 'zoh');

























function R = RotAxisAngle(n,phi) 
%
% OUTPUT:
%
% R: rotation matrix (active rotation)
%
% INPUT:
%
% n: rotation axis (also not normalized)
% phi: rotation angle (in radians)

% axis of rotation
n = n./norm(n);

N = [   0   -n(3)  n(2);
      n(3)     0  -n(1);
     -n(2)   n(1)    0 ];
 
 % rotation matrix
 R = eye(3) + sin(phi)*N + (1-cos(phi))*N*N;
 
end