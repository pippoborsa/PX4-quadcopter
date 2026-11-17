clear
close all
clc

syms t
syms phi real
syms phid real
syms phidd real
syms u real

L = 0.25; % [m]
Jxx = 0.0155; % [kg*m^2]
m = 1.07; % [kg]
c = 0.05; % [N/(m/s)]
Dz = 0.085; % [m] height of center of mass w.r.t. spherical hinge
g = 9.8; % [m/s^2]
k = 0.0183; % [N/PWM] force gain of each propeller

Jtot = Jxx + m*Dz^2; % transport formula for rotational inertia

T = (1/2)*(phid*Jtot*phid);                     % kinetic energy
% phi = 0 corresponds to vertical position
V = m*g*cos(phi)*Dz;                            % potential energy
D = (1/2)*(phid*c*phid);                        % dissipation function
Q = u*k*L;                                      % generalized force

% set free coordinates

q = [phi]';
qd = [phid]';
qdd = [phidd]';

%%%%=======================================================================
%%%%     Build Lagrange equations
%%%%     D(dT/dqd)/Dt - dT/dq + dV/dq - Q = 0 
%%%%=======================================================================

% Compute dT/dqd
dT_dqd = simplify(jacobian(T,qd))';

% Compute D(dT/dqd)/Dt
DdT_Dtdqd = simplify(jacobian(dT_dqd,t)) ...
            + simplify(jacobian(dT_dqd,q))  * qd ...
            + simplify(jacobian(dT_dqd,qd)) * qdd;

% note: the jacobian of T w.r.t. qd depends explicitly on time, q(t) and qd(t)

% Compute dT/dq
dT_dq = simplify(jacobian(T,q))' ;    
    
% Compute dV/dq
dV_dq = simplify(jacobian(V,q))' ; 

% Compute dD/dq
dD_dqd = simplify(jacobian(D,qd))' ; 
       
% Set up equations
Equations = simplify(DdT_Dtdqd - dT_dq + dD_dqd + dV_dq - Q)

% note: these equations refer to the CONTROLLED system!

% linearize the system and transform it to state space form

x1 = phi;   % x1 is position (phi)
x2 = phid;  % x2 is velocity (phid)

phidd_eq = solve(Equations, phidd);

% x1_dot = x2
% x2_dot = x2_ddot...

state_eqns = [x2; phidd_eq];

A = jacobian(state_eqns, [x1, x2]); % jacobian w.r.t. state
B = jacobian(state_eqns, u);        % jacobian w.r.t. control

A_lin_un = double(subs(A, [x1, x2, u], [0, 0, 0]));
B_lin_un = double(subs(B, [x1, x2, u], [0, 0, 0]));
C_lin_un = eye(2);
D_lin_un = 0;
disp('eigenvalues of uncontrolled linearized system')
eig(A_lin_un)

model_to_simulink = matlabFunction(state_eqns, 'file', 'InvertedPendulum','vars', [x1, x2, u]);
linear_sys = ss(A_lin_un, B_lin_un, C_lin_un, D_lin_un); % Create state-space model

Q = [1, 0; 0, 1];
R = 1;
N = [0; 0];

[K,S,P] = lqr(A_lin_un,B_lin_un,Q,R,N);
K


