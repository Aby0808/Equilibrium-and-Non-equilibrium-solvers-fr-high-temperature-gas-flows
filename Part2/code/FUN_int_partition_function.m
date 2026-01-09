function Q_int = FUN_int_partition_function(species, T, Species_data, Spectroscopic_data)
% Constants
k_B = 1.380649e-23; % Boltzmann constant (J/K)
h = 6.62607015e-34; % Planck constant (J s)
N = 3;              % number of terms for electronic states

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
r_e = Species_data.bond_length(species_idx1);  % bond length in meters
Tv = Species_data.frequency(species_idx1);  % charecteristic vibrational temperature
sigma = Species_data.symmetry_factor(species_idx1);  % symmetry number
g_elec = Spectroscopic_data.g_elec(elec_idx,1:N);  % Electronic degeneracy array
Te = Spectroscopic_data.energies(elec_idx,1:N);  % charecteristic electronic temperature

% For single atoms (like N, O, etc.), we don't have m2, r_e, or nu. Set them to 0 if needed.
if m2 == 0
    m2 = m1;  % For atoms, treat it as a single atom (not diatomic)
    r_e = 0;  % No bond length for atoms
    Tv = 0;   % No vibrational frequency for atoms
    sigma = 1;  % Symmetry number for single atoms
end

% Reduced mass for rotational partition function
mu = (m1 * m2) / (m1 + m2);  % reduced mass in kg

% Moment of inertia for diatomic molecules
I = mu * r_e^2;  % moment of inertia in kg m^2

% Rotational partition function
if r_e > 0
    % For diatomic molecules
    B = h / (8 * pi^2 * I);  % rotational constant in J
    Q_rot = (k_B * T) / (sigma * h * B);  % Rotational partition function
else
    % For single atoms, Q_rot = 1 (no rotational degrees of freedom)
    Q_rot = 1;
end

% Vibrational partition function (for diatomic molecules)
if Tv > 0
    Q_vib = 1 / (1 - exp(-Tv / T));  % Vibrational partition function
else
    Q_vib = 1;  % For single atoms (no vibrational degrees of freedom)
end

% Electronic partition function using spectroscopic data

% Te = Te * 1.98630e-23;  % cm^-1 to Joules conversion (1 cm^-1 = 1.98630e-23 J)
Q_elec = sum(g_elec .* exp(-Te / T));  % Sum over all electronic states

% Total internal partition function
Q_int = Q_rot * Q_vib * Q_elec;

end