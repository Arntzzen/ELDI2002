%% Block 6

%% Checking if our system is controlable and observable
%%%%%%%% Loading A- and B-matrices from Block 2 %%%%%%%%%%%%%%%%%%%%%%
A
B
Control = [B A*B]
Rank_Control = rank(Control)
det_Control = det(Control)

C = [0 1]
O = [C; C*A]
Rank_O = rank(O)
det_O = det(O)

 %% Designing full-state feedback controller
zeta = 0.591155033798898;
wn = 845.8018140975397;

syms k1 k2 lambda
k = [k1 k2];
system = A-B*k;
eq = charpoly(system, lambda)
sol = solve(eq, k);
[coef, powers] = coeffs(eq, lambda);
collected = arrayfun(@(c) collect(c, [k1 k2]), coef);
final = sum (collected .* lambda.^powers);
terms = collected .* lambda.^powers;
%arrayfun(@(t) pretty(t), terms);
pretty(terms(2))

first_term = collected(2)
second_term = collected(3)

two_zeta_wn = 2 * zeta * wn
wn_square = wn^2

eq1 = first_term == two_zeta_wn;
eq2 = second_term == wn_square;

sol_k = solve([eq1, eq2], [k1, k2]);
k1 = double(sol_k.k1)
k2 = double(sol_k.k2)

%{
k = [sol.k1, sol.k2]
eig(A-B*k)
double(eig(A-B*k))
%}