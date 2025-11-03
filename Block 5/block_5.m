%% Block 5

PO = 10;
Ts = 0.002;

zeta = (-log(PO/100)) / (sqrt(pi^2 + (log(PO/100))^2))
wn = 1 / (zeta * Ts)

gamma = 7788.62;
Kp = 0.0087;
Ki = gamma * Kp;