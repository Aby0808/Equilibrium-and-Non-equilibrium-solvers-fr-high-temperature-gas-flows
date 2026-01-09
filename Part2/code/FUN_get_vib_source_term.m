function [Omega] = FUN_get_vib_source_term(T, Tv, P, rho, y)

% this function calculates the vibrational energy source term for the vibrational energy equation
global MDOT SP_AVAIL

% calculating vibrational and electronic contribution to internal energy
% and Cvv
evib = FUN_get_mix_thermo_prop(SP_AVAIL,y,Tv,'e');
evib_eq = FUN_get_mix_thermo_prop(SP_AVAIL,y,T,'e');
Cvv = FUN_get_mix_thermo_prop(SP_AVAIL,y,Tv,'Cv');

Tau = FUN_get_tau(T,P,y);
Omega_v = rho*(evib_eq - evib)/Tau; % contribution from vibrational non equilibrium
w = FUN_kinetics_solver(T,Tv,rho,y);
Omega_sp = FUN_get_thermo_prop('N2',Tv,'e')*w(1) + ...
    FUN_get_thermo_prop('N',Tv,'e')*w(2);   % contribution from the production of molecular species in vibrational no eq

Omega = (Omega_v - Omega_sp)/(MDOT*Cvv);

end

function [Tau_tot] = FUN_get_tau(T,p,y)
% this function calculates the vibrational relaxation time

global K_B NA_AVO;

a=[221 180];
b=[0.029 0.0262];
m1 = FUN_get_thermo_prop('N2',T,'mol_wt');
m2 = FUN_get_thermo_prop('N',T,'mol_wt');
m=(m1*m2)/(m1+m2);
m=m/NA_AVO;
x = FUN_mass2mole_frac(y);

Tauv = exp(a.*(T^(-1/3) - b) - 18.421);
Tauv = sum(x.*Tauv)/(p/101325);

sigmav = 10^-21 * (50000/T)^2;

Tau_tot = Tauv + ((p/(K_B*T))*sqrt((8*K_B*T)/(pi*m))*sigmav)^-1;
end