function Kp = FUN_get_Kp_eq_stat(T)

global K_B NA_AVO SPECIES_DATA SPEC_DATA

P=101325;
V = 1;

eN = 470.8; %kJ/mol
eN2 = 0;

del_e1 = 1000*(2*eN - eN2)/NA_AVO;

%reaction 1 N2 <-> 2N
QN2 = FUN_total_partition_function('N2', T, V, SPECIES_DATA, SPEC_DATA);
QN = FUN_total_partition_function('N', T, V, SPECIES_DATA, SPEC_DATA);
Kp = exp(-del_e1/(K_B*T))*(QN^2 / QN2)*(K_B*T/P);

end