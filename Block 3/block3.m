
% Warm up task 1:

PO = 10;                                                % Percent overshoot
Tp = 0.005;                                             % Peak time

z = (-log(PO/100)) / (sqrt(pi^2 + (log(PO/100))^2));    % Calculated formula for zeta
wn = pi / (Tp*sqrt(1-z^2));                             % Calculated formula for wn

s = tf('s');
G = (wn^2) / (s^2 + 2*z*wn*s + wn^2);                   % Transfer function from warm up task 1

t = 0:0.00000001:0.02;                                  % Time vector
% y = step(G, t);                                       % Function generating step response for the transfer function "G" and time interval "t"

%{
% Plot the step response
figure;
plot(t, y);
xlabel('Time (s)');
ylabel('Response');
title('Step Response of the System');
grid on;

% Calculate the settling time and rise time
format long
settlingTime = stepinfo(G).SettlingTime
riseTime = stepinfo(G).RiseTime;
%}



%{
% Warm up task 2:

PO = 5;                                                      % Percent overshoot
Tp = 0.002;                                                  % Peak time
    
z = (-log(PO/100)) / (sqrt(pi^2 + (log(PO/100))^2));         % Calculated formula for zeta
wn = pi / (Tp*sqrt(1-z^2));                                      % Calculated formula for wn
balance = 10*z*wn

gam = 500;
%a1 = 5000;

%s = tf('s');
%G1 = ((gam1/a1)*(s+a1)*(wn^2)) / ((s^2 + 2*z*wn*s + wn^2)*(s+gam1));   % Transfer function from warm up task 1

t = 0:0.000001:0.005;                                             % Time vector
%y = step(G, t);                                                 % Function generating step response for the transfer function "G" and time interval "t"

a_val = [6500, 1000, 500];
gam_val = [8000, 1000, 500];




% Plot of the step response
figure;
hold on;
for i = 1:length(gam_val)
    %a = a_val(i);
    gam = gam_val(i);
    s = tf('s');
    G = ((gam/a)*(s+a)*(wn^2)) / ((s^2 + 2*z*wn*s + wn^2)*(s+gam));
    [y, ~] = step(G, t);
    plot(t, y, 'DisplayName', ['gamma = ' num2str(gam)]);
end

hold off;
legend show;
xlabel('Time (s)');
ylabel('Response');
title('Step Response of the System');
grid on;

% Calculate the settling time and rise time
format long
settlingTime = stepinfo(G).SettlingTime;
riseTime = stepinfo(G).RiseTime;
settlingValue = y(end);

%}


syms vC iL u R L C V0 G I0_nom a11 a12 a21 a22 b1 b2;

F1 = -R/L*iL - (1-u)/L*vC + V0/L;
F2 = -G/C*vC + (1-u)/C*iL - I0_nom/C;

Jac_AS =[a11 a12
         a21 a22];

Jac_BS = [b1
          b2];



% Transfer function of our system
syms IL VC U s;
X = [IL; VC];

%s = tf('s');
%eq = s*X == A*X + B*U
eq = s*X == Jac_AS*X + Jac_BS*U

%sol = solve(eq, [IL, VC])
sol = solve(eq, [IL VC]);

IL_col = collect(sol.IL, U);
disp('Transfer function:');
%pretty(IL_col);


%sim = vpa(simplify(sol.IL), 2)
%pretty(sim)

eq_val = s*X == A*X + B*U;
sol_eq_val = solve(eq_val, [IL, VC]);





% State-space system with transfer function
s = tf('s');

C_m = [1, 0];
D_m = 0;

sys_ss = ss(A, B, C_m, D_m)
sys_tf = tf(sys_ss)

%step(sys_tf);

pol = pole(sys_tf);
zer = zero(sys_tf);

K = evalfr(sys_tf, 0);
alpha = real(pol(1))

first_order = K * (1 / (s + abs(alpha)))

step(first_order);




% Closed-Loop
syms Kp Ki K aa s Wn_cl zeta_cl
G_hat = K * (1 / (s + aa));
Gc = Kp*((s+Ki/Kp)/s);

% Define the closed-loop transfer function
%tauc = Ki / Kp;
%Gcl = Gc*G_hat;
%Gcl = collect(Gcl, K);
%pretty(Gcl);

PI_cl = (G_hat * Gc) / (1 + G_hat * Gc);
PI_closed_loop = collect(PI_cl, s);


disp('PI Closed Loop:')
pretty(PI_closed_loop)

K = 3687.9;
aa = 301.3;

Ki_cl_sym = Wn_cl^2 / K
Kp_cl_sym = (2 * z * Wn_cl - aa) / K
s = -Ki_cl_sym / Kp_cl_sym;

sol = solve(abs((-Ki_cl_sym / Kp_cl_sym)) >= 10*z*Wn_cl, Wn_cl, 'ReturnConditions', true)

vpa(sol.conditions(1,1))
vpa(sol.conditions(2,1))
vpa(sol.conditions(3,1))


% Transfer function of the PI-controller
%alpha = 301.3;
wn_cl = 290;

Ki_cl = wn_cl^2 / K
Kp_cl = (2 * z * wn_cl - aa) / K
S = Ki_cl/Kp_cl;

sys_PI = tf([K*Kp_cl, Ki_cl*K],[1, (aa+K*Kp_cl), K*Ki_cl])
disp('sys_PI:')
stepinfo(sys_PI)
step(sys_PI)
