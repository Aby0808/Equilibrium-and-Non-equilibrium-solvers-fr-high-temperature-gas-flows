function [value] = FUN_get_thermo_prop(species, T, property)

%this function computes the thermodynamic property if the species
%species : name of the species
%property to be calculated

% Declare global constants
global RU NA_AVO SPECIES_DATA SPEC_DATA

% Parameters for the electronic states
Num = 3; % Number of terms for electronic states

% Find the index of the species in Species_data and Spectroscopic_data
species_idx1 = find(strcmp(SPECIES_DATA.species_names, species));

if isempty(species_idx1)
    error('Species %s not found in Species_data', species);
end

elec_idx = find(strcmp(SPEC_DATA.species_names, species));

if isempty(elec_idx)
    error('Species %s not found in Spectroscopic_data', species);
end

% Get species-specific data
m1 = SPECIES_DATA.mass(1, species_idx1);  % Mass of atom 1 (kg)
m2 = SPECIES_DATA.mass(2, species_idx1);  % Mass of atom 2 (kg) or 0 for single atom
Tv = SPECIES_DATA.frequency(species_idx1);  % Characteristic vibrational temperature
g_elec = SPEC_DATA.g_elec(elec_idx, 1:Num);  % Electronic degeneracy array
Te = SPEC_DATA.energies(elec_idx, 1:Num);  % Characteristic electronic temperature

% Molecular weight and specific gas constant
mol_wt = (m1 + m2) * NA_AVO; % Molecular weight in kg/mole
R_sp = RU / mol_wt;      % Specific gas constant (J/kg/K)

% Vibrational contribution to Cv
if m2==0
    Cv_vib=0;
else
    Cv_vib = ((Tv / T)^2 * exp(Tv / T) * R_sp) / (exp(Tv / T) - 1)^2;
end

% Electronic contribution to Cv
Z_elec = sum(g_elec .* exp(-Te / T)); % Partition function for electronic states
dlnZ_dT = (sum(g_elec .* Te .* exp(-Te / T)) / (T^2 * Z_elec)); % Temperature derivative
Cv_elec = RU * dlnZ_dT; % Contribution from electronic states

% Translational and rotational contributions to Cv
if m2==0
    Cv_trans_rot = 1.5 * R_sp;
else
    Cv_trans_rot = 2.5 * R_sp;
end

% Total Cv
Cv_total = Cv_trans_rot + Cv_vib + Cv_elec;

if strcmp(property, 'Cv')
    value = Cv_total;
elseif strcmp(property, 'Cp')
    % Cp = Cv + R_sp
    value = Cv_total + R_sp;
elseif strcmp(property, 'gamma')
    % Gamma = Cp / Cv
    Cp = Cv_total + R_sp;
    value = Cp / Cv_total;
elseif strcmp(property, 'R_sp')
    % Specific gas constant
    value = R_sp;
elseif strcmp(property, 'mol_wt')
    value = mol_wt;
elseif strcmp(property,'e')
    if m2==0
        value=0 + dlnZ_dT*R_sp*T^2;
    else
        value = ((Tv / T) / (exp(Tv / T) - 1))*R_sp*T  + dlnZ_dT*R_sp*T^2;
    end
else
    error('Invalid property specified: %s', property);
end
end
