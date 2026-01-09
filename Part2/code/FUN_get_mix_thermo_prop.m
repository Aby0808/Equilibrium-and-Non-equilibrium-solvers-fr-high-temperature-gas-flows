function [value] = FUN_get_mix_thermo_prop(species_array, mass_frac, T, property)

% this function computes the thermodynamic property of the mixture
%species_array: array containing the species
%prop: property to be found
%mass_frac: mass fraction of the mixture

nb_sp = size(species_array,2);

value=0;
if strcmp(property, 'Cv')
    for i=1:nb_sp
        Cv = FUN_get_thermo_prop(species_array(i),T,property);
        value = value + mass_frac(i)*Cv;
    end
elseif strcmp(property, 'Cp')
    for i=1:nb_sp
        Cp = FUN_get_thermo_prop(species_array(i),T,property);
        value = value + mass_frac(i)*Cp;
    end
elseif strcmp(property, 'gamma')
    Cp = 0;
    Cv = 0;
    for i=1:nb_sp
        Cv = Cv + mass_frac(i)*FUN_get_thermo_prop(species_array(i),T,'Cv');
        Cp = Cp + mass_frac(i)*FUN_get_thermo_prop(species_array(i),T,'Cp');
    end
    value = Cp/Cv;
elseif strcmp(property, 'R_sp')
    for i=1:nb_sp
        R = FUN_get_thermo_prop(species_array(i),T,property);
        value = value + mass_frac(i)*R;
    end
elseif strcmp(property, 'mol_wt')
    for i=1:nb_sp
        value = value + (mass_frac(i)/FUN_get_thermo_prop(species_array(i),T,property));
    end
    value = 1/value;
elseif strcmp(property, 'e')
    for i=1:nb_sp
        e = FUN_get_thermo_prop(species_array(i),T,property);
        value = value + mass_frac(i)*e;
    end
else
    error('Invalid property specified: %s', property);
end
end