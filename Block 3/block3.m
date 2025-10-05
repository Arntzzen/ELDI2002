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

gam = 5000;
%a1 = 5000;

%s = tf('s');
%G1 = ((gam1/a1)*(s+a1)*(wn^2)) / ((s^2 + 2*z*wn*s + wn^2)*(s+gam1));   % Transfer function from warm up task 1

t = 0:0.000001:0.01;                                             % Time vector
%y = step(G, t);                                                 % Function generating step response for the transfer function "G" and time interval "t"

    a_val = [5000, 1000, 500];



% Plot of the step response
figure;
hold on;
for i = 1:length(a_val)
    a = a_val(i);
    s = tf('s');
    G = ((gam/a)*(s+a)*(wn^2)) / ((s^2 + 2*z*wn*s + wn^2)*(s+gam));
    [y, ~] = step(G, t);
    plot(t, y, 'DisplayName', ['a = ' num2str(a)]);
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