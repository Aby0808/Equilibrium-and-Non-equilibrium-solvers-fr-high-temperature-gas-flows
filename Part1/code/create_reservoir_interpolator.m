function [enthalpy_interpolator, entropy_interpolator] = create_reservoir_interpolator(data)
    % Inputs:
    % - data: a table containing columns for T (Temperature), p (Pressure), Enthalpy, and Entropy.
    % Outputs:
    % - enthalpy_interpolator: a function handle for interpolating enthalpy based on T and p.
    % - entropy_interpolator: a function handle for interpolating entropy based on T and p.

    % Extract data columns
    T = data(:,1);           % Temperature column
    p = data(:,2);           % Pressure column
    h = data(:,5);    % Enthalpy column
    s = data(:,6);     % Entropy column

    % Create linear interpolators
    % Input points: [T, p]
    % Output values: h (for enthalpy) and s (for entropy)
    fprintf('\ncreating interpolator for enthalpy....')
    enthalpy_interpolator = scatteredInterpolant(T, p, h, 'linear', 'none'); % 'none' disables extrapolation
    fprintf('\ndone!')
    fprintf('\ncreating interpolator for entropy....')
    entropy_interpolator = scatteredInterpolant(T, p, s, 'linear', 'none');
    fprintf('\ndone!')
end
