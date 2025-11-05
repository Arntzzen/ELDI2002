%% Block 5

PO = 10;
Ts = 0.002;

zeta = (-log(PO/100)) / (sqrt(pi^2 + (log(PO/100))^2));
wn = 1 / (zeta * Ts);

gamma = 7788.62;
Kp = 0.0087;
Ki = gamma * Kp;


s1 = -499 + 682.9j;
a = 1549508692 / 158730;
co = 158730;
b = 602.5;
d = 126615937;

GclPI = co * ((s1 + a) / (s1 * ((s1)^2 + b*s1 + d)));


%%%%%%%%% Chat %%%%%%%%%%

% Finn total vinkel (i grader)
phi_total = angle(GclPI) * 180/pi;
%disp(phi_total)

phi_hoyre_side_av_likhetstegnet = -180 - phi_total;

hoyreside = tan(phi_hoyre_side_av_likhetstegnet)


syms gam
losning = 682.9 / (-499 + gam) == hoyreside;
l = double(solve(losning))