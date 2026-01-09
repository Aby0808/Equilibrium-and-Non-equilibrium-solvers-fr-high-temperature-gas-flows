function Q_tot = FUN_total_partition_function(species, T, V, Species_data, Spectroscopic_data)
% Constants
k_B = 1.380649e-23; % Boltzmann constant (J/K)
h = 6.62607015e-34; % Planck constant (J s)

% Find the index of the species in Species_data and Elec_spectral_data
species_idx1 = find(strcmp(Species_data.species_names, species));

if isempty(species_idx1)
    error('Species %s not found in Species_data', species);
end

% Get index of species in Elec_spectral_data
elec_idx = find(strcmp(Spectroscopic_data.species_names, species));

if isempty(elec_idx)
    error('Species %s not found in Elec_spectral_data', species);
end

% Get species-specific data
m1 = Species_data.mass(1,species_idx1);  % mass of atom 1 (kg)
m2 = Species_data.mass(2,species_idx1);  % mass of atom 2 (kg) or 0 for single atoms

% Reduced mass for translational partition function
m = m1 + m2;  % total mass in kg

Q_int = FUN_int_partition_function(species, T, Species_data, Spectroscopic_data);

Q_tr = ((2*pi*m*k_B*T)/h^2)^1.5 * V;

Q_tot = Q_tr*Q_int;

end