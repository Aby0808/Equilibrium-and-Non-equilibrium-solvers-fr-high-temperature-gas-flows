function [result] = process_nozzle_indirect_method(s0, h0, h_star, Area, A_x, A_star, ...
                                         F_rho_a_star, index_star, ...
                                         rbf_interpolator_rho, ...
                                         rbf_interpolator_speed, ...
                                         rbf_interpolator_pressure, ...
                                         rbf_interpolator_temperature)
% Routine to compute the inlet velocity using the bisection method,
% and calculate thermodynamic properties along the nozzle.
%
% Args:
% - s0: Initial entropy value.
% - h0: Reservoir enthalpy value.
% - Area: Table with area and position information.
% - A_x: Array of area values along the nozzle.
% - A_star: Throat area (minimum area).
% - F_rho_a_star: Mass flow rate divided by A_star (rho_star * a_star).
% - index_star: Index of the throat in A_x.
% - rbf_interpolator_*: Interpolators for thermodynamic properties.
%
% Returns:
% - Arrays of enthalpy, velocity, density, pressure, temperature, Mach number, and x positions.

% Extract x values from Area
% x_vec = Area.x(:);
x_vec = Area(:,1);

% Create two linspaces for enthalpy: from h0 to h_star and from h_star to 0
h_list_conv = linspace(h0, h_star, 1000);
h_list_div = linspace(h_star, 0, 1000)';

% Initialize result arrays
x_positions = [];
enthalpy_values = [];
velocity_values = [];
density_values = [];
pressure_values = [];
temperature_values = [];
mach_values = [];

% Split A_x into converging and diverging sections
A_x_conv = A_x(1:158);
A_x_div = A_x(159:end);

% Loop through converging section
for h = h_list_conv
    u = sqrt(2 * (h0 - h));
    input_interpolator = [log(s0), log(h)];
    rho = exp(rbf_interpolator_rho(input_interpolator));
    p = exp(rbf_interpolator_pressure(input_interpolator));
    T = exp(rbf_interpolator_temperature(input_interpolator));
    a = exp(rbf_interpolator_speed(input_interpolator));
    M = u / a;
    A = (F_rho_a_star / (rho * u)) * A_star;
    index = find_closest_index(A_x_conv, A);

    x = x_vec.x(index);

    x_positions = [x_positions; x];
    enthalpy_values = [enthalpy_values; h];
    velocity_values = [velocity_values; u];
    density_values = [density_values; rho];
    pressure_values = [pressure_values; p];
    temperature_values = [temperature_values; T];
    mach_values = [mach_values; M];
end

% Loop through diverging section
for h = h_list_div'
    u = sqrt(2 * (h0 - h));
    input_interpolator = [log(s0), log(h)];
    rho = exp(rbf_interpolator_rho(input_interpolator));
    p = exp(rbf_interpolator_pressure(input_interpolator));
    T = exp(rbf_interpolator_temperature(input_interpolator));
    a = exp(rbf_interpolator_speed(input_interpolator));
    M = u / a;
    A = (F_rho_a_star / (rho * u)) * A_star;
    index = find_closest_index(A_x_div, A);
    x = x_vec.x(index+158);

    if x < 0.05
        x_positions = [x_positions; x];
        enthalpy_values = [enthalpy_values; h];
        velocity_values = [velocity_values; u];
        density_values = [density_values; rho];
        pressure_values = [pressure_values; p];
        temperature_values = [temperature_values; T];
        mach_values = [mach_values; M];
    else
        break;
    end
end

result(:,1) = enthalpy_values;
result(:,2) = velocity_values;
result(:,3) = density_values;
result(:,4) = pressure_values;
result(:,5) = temperature_values;
result(:,6) = mach_values;
result(:,7) = x_positions;
end

function index = find_closest_index(array, value)
% Find the index of the closest value in the array to the given value
[~, index] = min(abs(array - value));
end
