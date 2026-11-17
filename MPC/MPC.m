clear
close all
clc

syms t
syms phi real
syms phid real
syms phidd real
syms u real

L = 0.25; % [m]
Jxx = 0.0159; % [kg*m^2]
m = 1.07; % [kg]
c = 0.011; % [N/(m/s)]
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

%%
Ts = 0.01;    % check se va bene online da review
dis_sys = c2d(linear_sys, Ts, 'zoh');

pr_horizon = 3;     % cambiare guardando simulazione circa 0.3/Ts
u_horizon = 1;

MV.min = -300;    % vincoli input tanto non dovrebbe sforare zona linearità
MV.max = 300;

OV(1).min = -pi/6;     
OV(1).max = pi/6;
OV(2).min = -pi;     
OV(2).max = pi;

W.ManipulatedVariables = 0.1/(MV.max)^2;
W.ManipulatedVariablesRate = 0.1*W.ManipulatedVariables;
W.OutputVariables = [50 3];



mpcobj = mpc(dis_sys, Ts, pr_horizon, u_horizon,W,MV,OV);


% Initial state setup
initialState = [-pi/6; 0];
Init_cond = mpcstate(mpcobj, initialState);
%%
% % Generate Explicit MPC
% range = generateExplicitRange(mpcobj);
% 
% % State limits (2 states: phi and phid)
% range.State.Min = [-pi; -10];
% range.State.Max = [pi; 10];
% 
% % Reference limits (2 outputs)
% range.Reference.Min = [-pi; -10];  % Adjust according to your specific application
% range.Reference.Max = [pi; 10];
% 
% % Manipulated variable limits (assuming 1 control input)
% range.ManipulatedVariable.Min = -200;
% range.ManipulatedVariable.Max = 200;
% 
% % Options for explicit MPC generation
% opt = generateExplicitOptions(mpcobj);
% opt.polyreduction = 100;
% 
% % Generate explicit MPC controller
% empcobj = generateExplicitMPC(mpcobj, range, opt);
