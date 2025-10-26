% Phase lead-compensator:
%clc; clear; close all;

% ---- plant + controller base values ----
K = 12.24;          % given plant gain
KI = 1230;         % integral gain used to set steady-state error (user value)
a = 301.3;          % plant pole location
% Target phase margin you want after compensation:
PM_target = 60;     % degrees
safety_deg = 5;     % safety margin (typical 3-10 deg)

% Define transfer functions
s = tf('s');
% Loop transfer L(s) (uncompensated) = K * KI / ( s*(s + a) )
L = (K*KI) / ( s*(s + a) );

% Plot uncompensated bode
figure; bode(L); grid on; title('Uncompensated loop L(j\omega)');


%% 1) find uncompensated gain crossover w_g (rad/s)
% Solve |L(jw)| = 1
syms w
eq = (K*KI)/(w*sqrt(w^2 + a^2)) == 1;   % magnitude of L(jw) == 1
w_g_sym = vpasolve(eq, w, 1);           % start guess near 1
w_g = double(w_g_sym);
w_g = w_g(w_g>0);

% If solver fails, fall back to numeric search
if isempty(w_g)
    fun = @(w) abs( (K*KI)/(1i*w*(1i*w + a)) ) - 1;
    w_g = fzero(fun, 1);
end
fprintf('Uncompensated gain crossover w_g = %.4f rad/s\n', w_g);

%% 2) uncompensated phase margin
% angle L(jw) = -90 - atan(w/a)  (in degrees)
phase_uncomp = -90 - atand(w_g/a);
PM_uncomp = 180 + phase_uncomp;
fprintf('Uncompensated phase at w_g = %.3f deg -> PM = %.3f deg\n', phase_uncomp, PM_uncomp);

%% 3) required additional phase lead (phi_m)
phi_needed = PM_target - PM_uncomp;    % how much more PM we need
phi_m = phi_needed + safety_deg;      % add safety margin
if phi_m <= 0
    warning('No additional phase needed (already >= target). phi_m <=0');
end
fprintf('Required extra phase (with safety) phi_m = %.3f deg\n', phi_m);

%% 4) compute alpha from sin(phi_m) = (alpha - 1)/(alpha + 1)
% use degrees trig functions
alpha = double( (1 + sind(phi_m)) / (1 - sind(phi_m)) );
% The dB boost the compensator gives at peak is 10*log10(alpha)
boost_dB = 10*log10(alpha);
fprintf('Computed alpha = %.4f -> compensator peak boost = %.3f dB\n', alpha, boost_dB);

%% 5) find frequency w_m where 20*log10(|L(jw)|) = -10*log10(alpha)
mag_target_dB = -boost_dB; % because compensator will add +boost_dB -> set crossover
mag_fun = @(w) 20*log10(abs(freqresp(L, w))) - mag_target_dB;

% freqresp(L,w) returns complex (MATLAB uses vector inputs; wrap)
% use fzero with a starting guess near w_g
try
    w_m = fzero(@(w) squeeze(20*log10(abs(freqresp(L, w)))) - mag_target_dB, w_g);
catch
    % if fzero fails, scan over logspace to find sign change
    ws = logspace(log10(max(1e-3,w_g/100)), log10(w_g*1000), 2000);
    vals = squeeze(20*log10(abs(freqresp(L, ws)))) - mag_target_dB;
    ix = find(vals(1:end-1).*vals(2:end) < 0, 1, 'first');
    if isempty(ix)
        error('Could not find w_m by scan. Try different starting guess or check L.');
    end
    w_m = (ws(ix)+ws(ix+1))/2;
end
fprintf('Chosen w_m = %.4f rad/s\n', w_m);

%% 6) compute T and uncompensated lead network
T = 1 / ( w_m * sqrt(alpha) );
% build uncompensated form of lead (note: include 1/alpha normalization)
Gc_base = (1 + alpha*T*s) / ( alpha * (1 + T*s) )   % no extra scalar Kc yet
figure; bode(Gc_base); grid on; title('Nominal lead G_c (before K-adjust)');
%{
%% 7) adjust compensator scalar so |Gc * L| = 1 at w_m
Ltemp = Gc_base * L;
mag_Ltemp_at_wm = abs(evalfr(Ltemp, 1i*w_m));
Kc = 1 / mag_Ltemp_at_wm;    % multiply Gc_base by this to make unity gain at w_m
Gc = Kc * Gc_base;
fprintf('Compensator scalar Kc = %.4f\n', Kc);

%% final compensated loop and margins
L_comp = Gc * L;

figure;
bode(L, L_comp, {1e-1, 1e5}); grid on;
legend('Uncompensated', 'Compensated');
title('Loop Bode: Uncompensated vs Compensated');

% show margins numerically
[Gm_u, Pm_u, Wcg_u, Wcp_u] = margin(L);
[~, Pm_c, ~, Wcp_c] = margin(L_comp);
fprintf('Uncompensated PM = %.3f deg at wc = %.4f rad/s\n', Pm_u, Wcp_u);
fprintf('Compensated PM   = %.3f deg at wc = %.4f rad/s\n', Pm_c, Wcp_c);

% Display final compensator
disp('Compensator Gc(s):');
Gc

% Optionally: plot closed-loop step to check time domain
Tcl = feedback(L_comp, 1);
figure; step(Tcl); grid on; title('Closed-loop step response with lead compensator');

stepinfo(Tcl)






% Phase lag-compensator:
%clc; clear; close all;

K = 12.24;
KI = 1230.8
a = 301.3;

PM_target = 60;
boost_dB = 10;

s = tf('s');
L = (K*KI) / ( s*(s + a) )

%% --- 1) Plot uncompensated loop ---
figure; bode(L); grid on; title('Uncompensated L(j\omega)');

%% --- 2) Gain crossover frequency (|L(jw)| = 1) ---
syms w
eq = (K*KI)/(w*sqrt(w^2+a^2)) == 1;
w_g = double(solve(eq,w));
w_g = w_g(w_g>0);
fprintf('Uncompensated gain crossover w_g = %.4f rad/s\n', w_g)


%% --- 3) Compute uncompensated phase margin ---
phase_uncomp = angle(evalfr(L, 1i*w_g)) * (180/pi)
PM_uncomp = 180 + phase_uncomp
fprintf('Uncompensated phase = %.3f deg, PM = %.3f deg\n', phase_uncomp, PM_uncomp)

%% --- 4) Choose alpha (low-frequency boost) ---
alpha = 10^(boost_dB/20)   % convert from dB to linear
fprintf('Chosen boost = %.2f dB  -->  alpha = %.4f\n', boost_dB, alpha)

%% --- 5) Choose lag pole/zero placement ---
% Place zero roughly one decade below crossover
w_z = w_g / 10         % zero frequency
T = 1 / w_z            % zero at 1/T
Gc = alpha * (1 + T*s) / (1 + alpha*T*s) %lag comp

%% --- 6) Form compensated open-loop transfer function ---
L_comp = Gc * L;

%% --- 7) Safe Bode plot limits ---
Wmin = max(1e-2, w_g/100);
Wmax = w_g*100;

figure;
bode(L, L_comp, {Wmin, Wmax});
grid on;
legend('Uncompensated', 'Compensated');
title('Bode Plot: Uncompensated vs Compensated');

%% --- 😎 Compute margins before and after ---
[Gm_u, Pm_u, Wcg_u, Wcp_u] = margin(L);
[Gm_c, Pm_c, Wcg_c, Wcp_c] = margin(L_comp);
fprintf('\nUncompensated: PM = %.2f° at ωc = %.4f rad/s\n', Pm_u, Wcp_u);
fprintf('Compensated:   PM = %.2f° at ωc = %.4f rad/s\n', Pm_c, Wcp_c);

%% --- 9) Closed-loop check ---
Tcl = feedback(L_comp, 1);
figure;
step(Tcl);
grid on;
title('Closed-loop step response with lag compensator');
stepinfo(Tcl)

z = 1/T;              % zero frequency (rad/s)
p = z/alpha;          % pole frequency (rad/s)
fprintf('Zero at ωz = %.3f rad/s, Pole at ωp = %.3f rad/s\n', z, p);
%}