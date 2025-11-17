%% Block 5
%%%%%%%%%%% Defining values and criterias%%%%%%%%%%%%
PO = 10;
Ts = 0.002;

zeta = (-log(PO/100)) / (sqrt(pi^2 + (log(PO/100))^2));
wn = 1 / (zeta * Ts)
Kv_comp = 1.7011e-7 * 5.8784e8
% These are values from our second order tf assigned to variables
s1 = -499 + 682.9j;
a = 1549508692 / 158730;
co = 158730;
b = 602.5;
d = 126615937;

% This is our loop tf with PI, but only with the s at origin and not
% the term  with gamma. We will calculate that later.
GolPI = co * ((s1 + a) / (s1 * ((s1)^2 + b*s1 + d)*(s1+20)))

%%%%%%%%% Angle criterion to find gamma %%%%%%%%%%
% Define desired pole sd
sd = -499+682.9j;

% calculating numeric values of each term (again, without gamma-term).
v1_num = sd + a;
v2_num = sd;
v3_num = sd + (-301.3+11264.4j);
v4_num = sd + (-301.3-11264.4j);
v5_num = sd + 20;

% Finding angles in degrees of each term.
v1_ang = angle(v1_num) * (180/pi)
v2_ang = angle(v2_num) * (180/pi)
v3_ang = angle(v3_num) * (180/pi)
v4_ang = angle(v4_num) * (180/pi)
v5_ang = angle(v5_num) * (180/pi)

% Calculating the sum of each angle and moving over to the right side
% of the equation and then getting rid of the arctan from the left side.
anglesum = v1_ang - v2_ang - v3_ang - v4_ang - v5_ang
rightside = -180 - anglesum
rightside_tan = tan(rightside * (pi/180))

% Solving for gamma.
syms gam
eq = 682.9 / (-499+gam) == rightside_tan
gamma = double(solve(eq))


%%%%%%%%% Magnitude criterion to find Kp and Ki %%%%%%%%%%

% Calculate the magnitude of the loop transfer function
magnitude = abs(GolPI)
% Set the desired magnitude for Kp
rightside = 1;
Kp = rightside / magnitude
Ki = gamma * Kp


%%%%%%%%%%%%%%% Plot root locus of Closed Loop PI %%%%%%%%%%%
desired_pole = [-499 + 682.9j, -499 - 682.9j]
s = tf('s');
Gol = ((Kp*(s+gamma)) / s) * (co*((s+a) / ((s^2 + b*s + d)*(s+20))))

ang_val = [v1_num, v2_num, v3_num, v4_num, v5_num];

% figure;
% rlocus(Gol);
% hold on;
% plot(desired_pole, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
%   plot(ang_val, 'b*')
% xlabel('Re(s)');
% ylabel('Im(s)');
% title('Root locus with desired poles');
% grid on;
% legend('Root Locus', 'Desired poles');



%%%%%%%%%%%%%%% Lag compensator %%%%%%%%%%%%%%%
al = (1 / (0.01*Kp*gamma*((co*a)/d)))% * 1.3977e-06
alpha = 1.7011e-7;
z = -499 / 50
p = z / alpha

comp = (s-z) / (s-p)
trans_func = comp * Gol
step(trans_func)

Ltf = trans_func / (1 + trans_func)
step(Ltf)
stepinfo(Ltf)

t = 0:0.01:1
ramp = t;               % Ramp (u(t) = t)



%%%%%%%%%%% Ramp input response plot %%%%%%%%%%
% --- Beregn systemrespons ---
[y, t_out] = lsim(Ltf, ramp, t);   % Lsim = lineær simulering

% --- Plot ---
figure;
plot(t_out, ramp, 'r--', 'LineWidth', 1.2); hold on;
plot(t_out, y, 'b', 'LineWidth', 1.5);
grid on;
xlabel('Tid (s)');
ylabel('Amplitude');
title('Ramp Input Response');
legend('Ramp input', 'System output');


%%%%%%%%%%%%%% Root Locus of Lag comp system %%%%%%%%%%%
figure;
rlocus(trans_func);
hold on;
plot(desired_pole, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
% plot(ang_val, 'b*')
xlabel('Re(s)');
ylabel('Im(s)');
title('Root locus with desired poles');
grid on;
legend('Root Locus', 'Desired poles');
% 
% % ess = double(limit(L, s, 0))