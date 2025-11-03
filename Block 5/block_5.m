%% Block 5

%% Defining variables
PO = 10;
Ts = 0.002;

zeta = (-log(PO/100)) / (sqrt(pi^2 + (log(PO/100))^2));
wn = 1 / (zeta * Ts);

gamma = -7788.62;       % Mulig feil fortegn?!?
Kp = 0.0087;
Ki = gamma * Kp;