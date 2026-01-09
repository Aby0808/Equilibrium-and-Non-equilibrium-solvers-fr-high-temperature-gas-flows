function [h, s] = get_reservoir_h_and_s(p, T, enthalpy_interpolator, entropy_interpolator)
    % This function returns the enthalpy and entropy for given pressure and temperature
    % using the provided linear interpolators.
    
    % Use the interpolators to get enthalpy and entropy at the specified p and T
    h = enthalpy_interpolator(T, p);
    s = entropy_interpolator(T, p);
end
