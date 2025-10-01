% Warm up task 1:

PO = 5;         % Percent overshoot
Tp = 0.002;     % Peak time

z = (-log(PO/100)) / (sqrt(pi^2 + (log(PO/100))^2));
wn = pi / (Tp*sqrt(1-z^2));

s = tf('s');
G = (wn^2) / (s^2 + 2*z*wn*s + wn^2);

