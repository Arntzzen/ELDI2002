V0 = 12;                                 % input voltage
Vc_nom = 11.9;                             % output voltage
D = (Vc_nom-V0)/Vc_nom;                  % Nominal duty cycle
P_nom = 100;                             % output power
I0_nom = P_nom/Vc_nom;                   % Output current
fs = 50000;                              % switching frequency
I0_max = 0.1 * I0_nom;                   % rippel current 10%
Vc_max = 0.05 * Vc_nom;                  % rippel voltage 5%
L_min = (V0*D)/(fs*I0_max);              % Min inductan
C_min = (I0_nom*D)/(fs*Vc_max);          % Min capacitance
L = 1.05 * L_min;
C = 1.05 * C_min;
n = 0.95;                                % Nominal efficiency

Il_nom = P_nom/(n*V0);                   % nominal input current

P_loss = (1-n)*P_nom;                    % Power loss

R = 3.75/Il_nom^2;                       % 0.75 * P_loss - Resistance
G = 1.25/Vc_nom^2;                       % 25% * P_loss - Conductance

Ilbar_ikkebruk = (V0 + sqrt(V0^2 - 4*R*(G*Vc_nom^2 + I0_nom*Vc_nom))) / (2*R); % Write in report: We remove this equation because of improbable values.
Ilbar = (V0 - sqrt(V0^2 - 4*R*(G*Vc_nom^2 + I0_nom*Vc_nom))) / (2*R);

u_ikkebruk = 1 - (G*Vc_nom + I0_nom) / Ilbar_ikkebruk; % Write in report: We remove this equation from the solution set (result = 0.9734) because of off-nominal value
u = 1 - (G*Vc_nom + I0_nom) / Ilbar;

Kp = 0.5;
KI = 0.5;

Kp_c = 500;
KI_c = 200;
%{
dILdt = -R/L*IL - (1-u)/L*Vc + V0/L
dVcdt = -G/C*Vc + (1-u)/C*IL - I0_nom/C
%}