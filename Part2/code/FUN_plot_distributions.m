function FUN_plot_distributions(species, x_array, T_array, Tv_array, x_query)
    % Function to plot rotational, vibrational, and electronic distributions in eV.
    % Separate plots are created for each distribution, showing all x-query locations.

    % Global variables
    global SPECIES_DATA SPEC_DATA

    % Ensure species exists in the dataset
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
    m2 = SPECIES_DATA.mass(2, species_idx1);  % Mass of atom 2 (kg)
    Theta_r = SPECIES_DATA.frequency(species_idx1);  % Rotational characteristic temp (K)
    Theta_v = SPECIES_DATA.frequency(species_idx1);  % Vibrational characteristic temp (K)
    g_elec = SPEC_DATA.g_elec(elec_idx, :);  % Electronic degeneracy array
    Theta_e = SPEC_DATA.energies(elec_idx, :);  % Electronic characteristic temperatures (K)

    % Interpolate T and Tv at specified x-locations
    T_query = interp1(x_array, T_array, x_query, 'linear', 'extrap');
    Tv_query = interp1(x_array, Tv_array, x_query, 'linear', 'extrap');

    % Constants
    k_B = 8.6173e-5; % Boltzmann constant in eV/K
    h = 4.1357e-15; % Planck's constant in eV·s
    c = 3e8; % Speed of light in m/s

    % Rotational quantum numbers
    J = 0:30; % Rotational quantum numbers
    g_J = 2 * J + 1; % Rotational degeneracy
    E_J = Theta_r * J .* (J + 1) * k_B; % Rotational energy in eV

    % Vibrational quantum numbers
    v = 0:10; % Vibrational quantum numbers
    E_v = Theta_v * (v + 0.5) * k_B; % Vibrational energy in eV

    % Electronic levels
    E_e = Theta_e * k_B; % Electronic energy levels in eV

    % Initialize data storage for plots
    rotational_distributions = [];
    vibrational_distributions = [];
    electronic_distributions = [];

    % Compute distributions for each x-location
    for i = 1:length(x_query)
        % Current interpolated temperatures
        T_loc = T_query(i); % Translational temperature
        T_v_loc = Tv_query(i); % Vibrational temperature

        % Rotational distribution
        f_r = g_J .* exp(-E_J / (k_B * T_loc));
        f_r = f_r / sum(f_r); % Normalize
        rotational_distributions = [rotational_distributions; f_r ./ g_J];

        % Vibrational distribution
        f_v = exp(-E_v / (k_B * T_v_loc));
        f_v = f_v / sum(f_v); % Normalize
        vibrational_distributions = [vibrational_distributions; f_v];

        % Electronic distribution
        f_e = g_elec .* exp(-E_e / (k_B * T_loc)); % Boltzmann factor using Theta_e
        f_e = f_e / sum(f_e); % Normalize
        electronic_distributions = [electronic_distributions; f_e ./ g_elec];
    end

    % Plot Rotational Distribution
    figure;
    hold on;
    for i = 1:length(x_query)
        plot(E_J, rotational_distributions(i, :), 'DisplayName', ['x = ', num2str(x_query(i)), ' m']);
    end
    xlabel('Rotational Energy (eV)');
    ylabel('X_J / g_J');
    grid on
    title(['Rotational Distribution ',species]);
    legend('show');
    hold off;

    % Plot Vibrational Distribution
    figure;
    hold on;
    for i = 1:length(x_query)
        plot(E_v, vibrational_distributions(i, :), 'DisplayName', ['x = ', num2str(x_query(i)), ' m']);
    end
    xlabel('Vibrational Energy (eV)');
    ylabel('X_v');
    grid on
    title(['Vibrational Distribution ',species]);
    legend('show');
    hold off;

    % Plot Electronic Distribution
    figure;
    hold on;
    for i = 1:length(x_query)
        plot(E_e, electronic_distributions(i, :), 'DisplayName', ['x = ', num2str(x_query(i)), ' m']);
    end
    grid on
    xlabel('Electronic Energy (eV)');
    ylabel('X_i / g_i');
    title(['Electronic Distribution ',species]);
    legend('show');
    hold off;
end
