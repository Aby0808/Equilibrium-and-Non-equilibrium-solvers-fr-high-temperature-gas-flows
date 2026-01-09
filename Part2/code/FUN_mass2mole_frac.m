function [x]=FUN_mass2mole_frac(y)
% this function converts mass fractions to mole fractions

global SP_AVAIL
T=0; % garbage value for T not required for conversion
Mmix = FUN_get_mix_thermo_prop(SP_AVAIL,y,T,'mol_wt');
% x = y .* ([Mmix/FUN_get_thermo_prop(SP_AVAIL(1),T,'mol_wt') Mmix/FUN_get_thermo_prop(SP_AVAIL(2),T,'mol_wt')])';
x(1) = y(1)*Mmix/FUN_get_thermo_prop(SP_AVAIL(1),T,'mol_wt');
x(2) = y(2)*Mmix/FUN_get_thermo_prop(SP_AVAIL(2),T,'mol_wt');
end