function [rbf_rho, rbf_speed, rbf_pressure, rbf_temperature] = construct_rbf_interpolators(Enthalpy, Entropy, rho, speed_of_sound, Pressure, Temperature)
    % Constructs custom RBF interpolators using Thin Plate Spline or other kernels.
    %
    % Inputs:
    %   Enthalpy, Entropy: Independent variables
    %   rho, speed_of_sound, Pressure, Temperature: Properties to interpolate
    %
    % Outputs:
    %   rbf_rho, rbf_speed, rbf_pressure, rbf_temperature: Interpolator functions

    % Prepare the input data
    x = log([Entropy, Enthalpy]); % Input: log of Entropy and Enthalpy

    % Kernel type and parameters
    kernelType = 'multiquadric'; % Options: 'thin_plate', 'cubic', 'quintic', etc.
    epsilon = 1.0; % Kernel scaling parameter
    regularization = 1e-6; % Regularization parameter

    % Fit RBF interpolators
    fprintf('\nCreating RBF interpolator for density....\n');
    rbf_rho = fitCustomRBF(x, log(rho), kernelType, epsilon, regularization);
    fprintf('..Done!\n');

    fprintf('\nCreating RBF interpolator for speed of sound....\n');
    rbf_speed = fitCustomRBF(x, log(speed_of_sound), kernelType, epsilon, regularization);
    fprintf('..Done!\n');

    fprintf('\nCreating RBF interpolator for pressure....\n');
    rbf_pressure = fitCustomRBF(x, log(Pressure), kernelType, epsilon, regularization);
    fprintf('..Done!\n');

    fprintf('\nCreating RBF interpolator for temperature....\n');
    rbf_temperature = fitCustomRBF(x, log(Temperature), kernelType, epsilon, regularization);
    fprintf('..Done!\n');
end

% Helper function to fit an RBF interpolator
function rbf_model = fitCustomRBF(x, y, kernelType, epsilon, regularization)
    % Custom RBF interpolator fitting
    %
    % Inputs:
    %   x: Nx2 matrix of independent variables
    %   y: Nx1 vector of dependent variable
    %   kernelType: RBF kernel ('thin_plate', 'cubic', etc.)
    %   epsilon: Scaling parameter for the kernel
    %   regularization: Regularization parameter
    %
    % Output:
    %   rbf_model: Function handle for interpolation

    % Compute pairwise distances between training points (manual implementation)
    n = size(x, 1);
    distances = sqrt(sum((reshape(x, n, 1, size(x, 2)) - reshape(x, 1, n, size(x, 2))).^2, 3));

    % Generate the kernel matrix
    phi = computeRBFKernel(distances, kernelType, epsilon);

    % Add regularization to the diagonal
    phi = phi + regularization * eye(n);

    % Solve for weights (coefficients)
    weights = phi \ y;

    % Return the RBF interpolator as a function handle
    rbf_model = @(xi) evaluateRBF(xi, x, weights, kernelType, epsilon);
end


% Function to compute the RBF kernel
function phi = computeRBFKernel(distances, kernelType, epsilon)
    % Compute the kernel matrix for the specified RBF type
    %
    % Inputs:
    %   distances: NxN pairwise distance matrix
    %   kernelType: RBF kernel type ('thin_plate', 'cubic', etc.)
    %   epsilon: Scaling parameter
    %
    % Output:
    %   phi: NxN kernel matrix

    switch lower(kernelType)
        case 'thin_plate'
            phi = distances.^2 .* log(distances + eps); % Thin Plate Spline
        case 'cubic'
            phi = distances.^3; % Cubic kernel
        case 'quintic'
            phi = -distances.^5; % Quintic kernel
        case 'gaussian'
            phi = exp(-(epsilon * distances).^2); % Gaussian kernel
        case 'multiquadric'
            phi = -sqrt(1 + (epsilon * distances).^2); % Multiquadric kernel
        case 'inverse_multiquadric'
            phi = 1 ./ sqrt(1 + (epsilon * distances).^2); % Inverse Multiquadric kernel
        otherwise
            error('Unsupported kernel type: %s', kernelType);
    end
end

% Function to evaluate the RBF at new points
function result = evaluateRBF(xi, x, weights, kernelType, epsilon)
    % Evaluate the RBF interpolator at new points
    %
    % Inputs:
    %   xi: Mx2 matrix of new points
    %   x: Nx2 matrix of training points
    %   weights: Nx1 vector of RBF weights
    %   kernelType: RBF kernel type
    %   epsilon: Scaling parameter
    %
    % Output:
    %   result: Mx1 vector of interpolated values

    % Compute pairwise distances between xi and x
    m = size(xi, 1);
    n = size(x, 1);
    distances = sqrt(sum((reshape(xi, m, 1, size(xi, 2)) - reshape(x, 1, n, size(x, 2))).^2, 3));

    % Compute the kernel values for the new points
    phi = computeRBFKernel(distances, kernelType, epsilon);

    % Compute the interpolated result
    result = phi * weights;
end
