function dydx = ode_system(x, prim_var)

global RU SP_AVAIL MDOT

% Extract variables
num_species = length(SP_AVAIL);
y = prim_var(1:num_species);
u = prim_var(num_species + 1);
T = prim_var(num_species + 2);
Tv = prim_var(num_species + 3);

% Mixture properties
R_mix = FUN_get_mix_thermo_prop(SP_AVAIL, y, T, 'R_sp');  %T is not involved in calculation of R
Cpmix = FUN_get_mix_thermo_prop(SP_AVAIL, y, sqrt(T*Tv), 'Cp');

% Compute derivatives
dydx = zeros(num_species + 2, 1);

detA = u^2 * Cpmix - R_mix*(T*Cpmix + u^2); % computing determinants

%get species production terms
rho = MDOT/u;
P = rho*R_mix*T;
omega_kinetic = FUN_kinetics_solver(T,Tv,rho,y);

%computing common terms
omi_Mi = 0;
hi_omi = 0;
for i = 1:num_species
    Cpi = FUN_get_thermo_prop(SP_AVAIL(i),sqrt(T*Tv),'Cp');
    mwi = FUN_get_thermo_prop(SP_AVAIL(i),T,'mol_wt');  %T is not involved in computation of mol wt
    omi_Mi = omi_Mi + (omega_kinetic(i)/mwi);
    hi_omi = hi_omi + Cpi*T*omega_kinetic(i);  % #### have to check T this #############333
end

% Species equations
for i = 1:num_species
    dydx(i) = omega_kinetic(i) / MDOT;
end

% Momentum equation
dydx(num_species + 1) = (-(u*RU*T*Cpmix*omi_Mi) + (u*R_mix*hi_omi))/(detA*MDOT);

% Energy equation
dydx(num_species + 2) = (-(u^2 - (R_mix*T))*hi_omi + u^2 * RU*T*omi_Mi)/(detA*MDOT);

% Vibrational energy equation
omega_vib = FUN_get_vib_source_term(T,Tv,P,rho,y);
dydx(num_species + 3) = omega_vib;

end