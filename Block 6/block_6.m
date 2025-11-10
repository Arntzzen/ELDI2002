%% Block 6

%%%%%%%% Checking if our system is controlable and observable %%%%%%%%
%%%%%%%% Loading A- and B-matrices from Block 2 %%%%%%%%%%%%%%%%%%%%%%
C = [B A*B]
Rank_C = rank(C)
det_C = det(C)

O = [C; C*A]
Rank_O = rank(O)
det_O = det(O)