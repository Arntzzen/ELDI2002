%{

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
settlingTime = stepinfo(G).SettlingTime
riseTime = stepinfo(G).RiseTime
settlingValue = y(end)

%}


syms vC iL u R L C V0 G I0_nom a11 a12 a21 a22 b1 b2;

F1 = -R/L*iL - (1-u)/L*vC + V0/L;
F2 = -G/C*vC + (1-u)/C*iL - I0_nom/C;

dF1diL = diff(F1, iL);
dF1dvC = diff(F1, vC);
dF2diL = diff(F2, iL);
dF2dvC = diff(F2, vC);

Jac_A = [dF1diL dF1dvC
         dF2diL dF2dvC];

Jac_AS =[a11 a12
         a21 a22]



dF1du = diff(F1, u);
dF2du = diff(F2, u);

Jac_B = [dF1du
         dF2du];
Jac_BS = [b1
          b2]



% Transfer function of our system
syms IL VC U s;
X = [IL; VC]

%s = tf('s');
%eq = s*X == A*X + B*U
eq2 = s*X == Jac_AS*X + Jac_BS*U

%sol = solve(eq, [IL, VC])
sol2 = solve(eq2, [IL VC])

IL_expr = sol2.IL

disp('Y(s):');
pretty(IL_expr);

IL_col = collect(IL_expr, U);
disp('Transfer function:')
pretty(IL_col)

latex(sol2.IL);


%sim = vpa(simplify(sol.IL), 2)
%pretty(sim)

%pretty(sol2.IL)