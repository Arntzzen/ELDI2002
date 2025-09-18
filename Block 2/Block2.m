%Defining variables
V0 = 12;                                 % input voltage
Vc_nom = 16;                             % output voltage
D = (Vc_nom-V0)/Vc_nom;                  % Nominal duty cycle
P_nom = 100;                             % output power
I0_nom = P_nom/Vc_nom;                   % Output current
fs = 50000;                              % switching frequency
I0_max = 0.1 * I0_nom;                   % rippel current 10%
Vc_max = 0.05 * Vc_nom;                  % rippel voltage 5%
L_min = (V0*D)/(fs*I0_max);              % Min inductan
C_min = (I0_nom*D)/(fs*Vc_max);          % Min capacitance
L = 1.05 * L_min;                        % Chosen value for L
C = 1.05 * C_min;                        % Chosen value for C
n = 0.95;                                % Nominal efficiency
Il_nom = P_nom/(n*V0);                   % nominal input current
P_loss = (1-n)*P_nom;                    % Power loss
R = 3.75/Il_nom^2;                       % 0.75 * P_loss - Resistance
G = 1.25/Vc_nom^2;                       % 25% * P_loss - Conductance


%Equilibrium values
Ilbar = (V0 - sqrt(V0^2 - 4*R*(G*Vc_nom^2 + I0_nom*Vc_nom))) / (2*R);
u = 1 - (G*Vc_nom + I0_nom) / Ilbar;


% Inner-loop PI-values
Kp = 0.5;
KI = 0.5;


% Outer-loop PI-values
Kp_c = 500;
KI_c = 200;


% Differential equations of our system
% dILdt = -R/L*IL - (1-u)/L*Vc + V0/L           %IL is non-defined variable and therefore commented out
% dVcdt = -G/C*Vc + (1-u)/C*IL - I0_nom/C       %Vc is non-defined variable and therefore commented out


% Jacobian for A-matrix
syms vC iL u;

F1 = -R/L*iL - (1-u)/L*vC + V0/L;
F2 = -G/C*vC + (1-u)/C*iL - I0_nom/C;

dF1diL = diff(F1, iL);
dF1dvC = diff(F1, vC);
dF2diL = diff(F2, iL);
dF2dvC = diff(F2, vC);

Jac_A = [dF1diL dF1dvC
         dF2diL dF2dvC]

% Evaluate the Jacobian matrix at the equilibrium point
Jac_A_eq = subs(Jac_A, {vC, iL, u}, {Vc_nom, Ilbar, u})


% Jacobian for B-matrix
syms vC iL u;

dF1du = diff(F1, u);
dF2du = diff(F2, u);

Jac_B = [dF1du
         dF2du]