PO = 10;
zeta = (-log(PO/100)) / (sqrt(pi^2 + (log(PO/100))^2));

alpha = -301.3;
K = 12.24;
wn_cl = 250;

Ki_cl = wn_cl^2 / K
Kp_cl = (2 * z * wn_cl - alpha) / K
S = Ki_cl/Kp_cl

sys_PI = tf([K*Kp_cl, Ki_cl*K,],[1, (alpha+K*Kp_cl), K*Ki_cl])
disp('sysbab:')
stepinfo(sys_PI)
step(sys_PI)

%syms wn a z s
%sys_symbolic = ((wn^2/a)*(s+a)*(wn^2) / ((s^2 + 2*z*wn*s) + wn^2))
%pretty(sys_symbolic)