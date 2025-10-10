PO = 10;
zeta = (-log(PO/100)) / (sqrt(pi^2 + (log(PO/100))^2));

alpha = 301.3;
K = 12.24;
wn_cl = 255;

Ki_cl = wn_cl^2 / K
Kp_cl = (2 * z * wn_cl - alpha) / K
S = Ki_cl/Kp_cl

s = tf('s');
sys_pu = ((1/alpha)*(s+alpha)*(wn_cl^2)) / ((s^2 + 2*zeta*wn_cl*s + wn_cl^2));
%step(sys_pu);
stepinfo(sys_pu)

sysbab = tf([K*Kp_cl, Ki_cl*K,],[1, (alpha+K*Kp_cl), K*Ki_cl])
disp('sysbab:')
stepinfo(sysbab)
step(sysbab)

%syms wn a z s
%sys_symbolic = ((wn^2/a)*(s+a)*(wn^2) / ((s^2 + 2*z*wn*s) + wn^2))
%pretty(sys_symbolic)