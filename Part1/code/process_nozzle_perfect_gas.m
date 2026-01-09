function [results_fr] = process_nozzle_perfect_gas(gamma, R, p0, T0, Area, A_x, A_star, index_star)
% Routine to compute flow properties in a nozzle for a perfect gas using the area-Mach number relation.
%
% Args:
% - gamma: Specific heat ratio (Cp/Cv).
% - R: Gas constant (J/kg·K).
% - p0: Total (stagnation) pressure (Pa).
% - T0: Total (stagnation) temperature (K).
% - Area: Array or table with area and position information (must have 'x' and 'A' columns).
% - A_x: List or array of area values along the nozzle.
% - A_star: Throat area (minimum area).
% - index_star: Index of the throat in A_x.
%
% Returns:
% - Lists of enthalpy, velocity, density, pressure, temperature, Mach number, and x positions.

    % Extract x and area data
    x_vec = table2array(Area(:, 1));  % Assuming Area is a table with x in the first column
    A_vec = table2array(Area(:, 2));  % Assuming Area has area data in the second column

    % Initialize result arrays
    enthalpy_values = [];
    velocity_values = [];
    density_values = [];
    pressure_values = [];
    temperature_values = [];
    mach_values = [];
    x_positions = [];

    % Calculate the area-to-throat area ratio for each position
    AAstar = A_vec / A_star;

    % Call the function to calculate the subsonic and supersonic Mach numbers
    [Msub, Msup] = FUN_Area_Mach_solve(AAstar, gamma);

    % Loop over the nozzle sections (converging and diverging)
    for i = 1:length(A_x)-1
        
        % Choose the appropriate Mach number (subsonic or supersonic)
        if i < index_star
            M = Msub(i); % Subsonic region
        elseif i == index_star
            M = 1;
        elseif i>index_star
            M = Msup(i+1); % Supersonic region
        end
        
        % Calculate density and temperature from the stagnation values and Mach number
        T = T0 / (1 + (gamma - 1) / 2 * M^2);
        p = p0 * (1 + (gamma - 1) / 2 * M^2) ^ (-gamma / (gamma - 1));
        
        % Using the ideal gas law: p = rho * R * T, calculate density
        rho = p / (R * T);
        
        % Calculate velocity (using the Mach number)
        u = M * sqrt(gamma * R * T);
        
        % Calculate enthalpy using the perfect gas relation: h = Cp * T
        h = gamma * R * T / (gamma - 1); % For perfect gas, Cp = gamma * R / (gamma - 1)
        
        % Append results
        x_positions = [x_positions; x_vec(i)];
        enthalpy_values = [enthalpy_values; h];
        velocity_values = [velocity_values; u];
        density_values = [density_values; rho];
        pressure_values = [pressure_values; p];
        temperature_values = [temperature_values; T];
        mach_values = [mach_values; M];
    end

    % Store results in a structure
    results_fr.enthalpy = enthalpy_values;
    results_fr.velocity = velocity_values;
    results_fr.density = density_values;
    results_fr.pressure = pressure_values;
    results_fr.temperature = temperature_values;
    results_fr.mach = mach_values;
    results_fr.x = x_positions;
end


function [Msub, Msup] = FUN_Area_Mach_solve(AAstar, g)
% Function to solve for subsonic and supersonic Mach numbers using the area-Mach relation.
% AAstar: Area ratio (A/A*) for each position along the nozzle.

    gm1 = g - 1;
    gp1 = g + 1;

    Msub = zeros(size(AAstar));  % Initialize subsonic Mach number array
    Msup = zeros(size(AAstar));  % Initialize supersonic Mach number array

    for i = 1:length(AAstar)-1
        % Define the area-Mach relation function
        area_mach_function = @(M) (1 / M^2) * ((2 + gm1 * M^2) / gp1)^(gp1 / gm1) - AAstar(i)^2;

        % Subsonic Mach number solution (expected to be between 0 and 1)
        try
            Msub(i) = fzero(area_mach_function, [1e-6, 1]);  % Subsonic Mach number
        catch
            disp(['Subsonic fzero failed for i = ', num2str(i), ', AAstar = ', num2str(AAstar(i))]);
            Msub(i) = NaN; % Assign NaN to flag the issue
        end

        % Supersonic Mach number solution (expected to be >1)
        try
            Msup(i) = fzero(area_mach_function, [1 + 1e-6, 50]);  % Supersonic Mach number
        catch
            disp(['Supersonic fzero failed for i = ', num2str(i), ', AAstar = ', num2str(AAstar(i))]);
            Msup(i) = NaN; % Assign NaN to flag the issue
        end
    end
end
