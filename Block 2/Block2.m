%Defining variables
V0 = 12;                                 % input voltage
vC_nom = 16;                             % output voltage
D = (vC_nom-V0)/vC_nom;                  % Nominal duty cycle
P_nom = 100;                             % output power
I0_nom = P_nom/vC_nom;                   % Output current
fs = 50000;                              % switching frequency
I0_max = 0.1 * I0_nom;                   % rippel current 10%
Vc_max = 0.05 * vC_nom;                  % rippel voltage 5%
L_min = (V0*D)/(fs*I0_max);              % Min inductan
C_min = (I0_nom*D)/(fs*Vc_max);          % Min capacitance
L = 1.05 * L_min;                        % Chosen value for L
C = 1.05 * C_min;                        % Chosen value for C
n = 0.95;                                % Nominal efficiency
iL_nom = P_nom/(n*V0);                   % nominal input current
P_loss = (1-n)*P_nom;                    % Power loss
R = 3.75/iL_nom^2;                       % 0.75 * P_loss - Resistance
G = 1.25/vC_nom^2;                       % 25% * P_loss - Conductance


%Equilibrium values
iL_bar = (V0 - sqrt(V0^2 - 4*R*(G*vC_nom^2 + I0_nom*vC_nom))) / (2*R);
u_bar = 1 - (G*vC_nom + I0_nom) / iL_bar;


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
         dF2diL dF2dvC];

% Evaluate the Jacobian matrix at the equilibrium point
Jac_A_lin = double(subs(Jac_A, {vC, iL, u}, {vC_nom, iL_bar, u_bar}));

A = Jac_A_lin;


% Jacobian for B-matrix
syms vC iL u;

dF1du = diff(F1, u);
dF2du = diff(F2, u);

Jac_B = [dF1du
         dF2du];

Jac_B_lin = double(subs(Jac_B, {vC, iL, u}, {vC_nom, iL_bar, u_bar}));

B = Jac_B_lin_sim;


% Calculations for x-matrix, x_bar-matrix and E-matrix
x = [iL
     vC];

x_ref = [iL_bar
         vC_nom];

E1 = -A * x_ref - B * u_bar;
E = double(E1);



% PI - Inner Loop
syms x_c KpIL KiIL;
e = [1; 0];
e_t = transpose(e);

x_3 = [x; x_c];
A_3 = [A-B*KpIL*e_t   B*KiIL
        -e_t         0 ]

B_3 = [B*KpIL
        1  ];

E_3 = [E
       0];

dx_3dt = A_3 * x_3 + B_3 * e_t*x_ref + E_3;



% PI - Outer Loop
syms x_c2 KpOL KiOL;
e2 = [0;1;0];
e2_t = transpose(e2);
x_refIL = e_t * x_ref;

x_4 = [x_3; x_c2];
A_4 = [A_3-B_3*KpOL*e2_t, B_3*KiOL; -e2_t, 0];

B_4 = [B_3*KpOL; 1];

E_4 = [E_3; 0];

dx_4dt = A_4 * x_4 + B_4 * e2_t*x_refIL + E_4



KiIL = 0.5



% Define the range of Kp -- > parameter to be swept
Kp_values = linspace(1, 5, 100); % for example, 10 steps from 1 to 10
% Prepare figure
figure;
hold on;
colormap(jet); % Color map for Kp
cmap = colormap;
nColors = size(cmap, 1);
% Loop over Kp values
for i = 1:length(Kp_values)
           KpIL = Kp_values(i);
          A_for = [A-B*KpIL*e_t   B*KiIL
        -e_t         0 ] % Example matrix, You put your own A matrix.
          % Compute eigenvalues
          eigvals = eig(A_for);
          % Choose color based on Kp
          colorIdx = round((i-1)/(length(Kp_values)-1)*(nColors-1)) + 1;
          plot(real(eigvals), imag(eigvals), 'x', 'Color', cmap(colorIdx,:), 'LineWidth', 1.5);
end
% Plot formatting
xlabel('Real Axis [s^{-1}]');
ylabel('Imaginary Axis [s^{-1}]');
title('Eigenvalue Sweep over Kp');
colorbar('Ticks', linspace(0, 1, 5), 'TickLabels', round(linspace(1, 5, 5), 1));
grid on;
axis equal;
%}