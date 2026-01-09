function [h_star, s_star] = compute_hstar_sstar(s0, h0, sos_interp)

% % this function calculates the s* and h* (conditions at throat), using the
% % relation sqrt(2h0 - h*) = a(h*,s0)
% % s0    : entropy at reservoir
% % h0    : enthaly at reservoir
% % sos_interp : interpolator of speed of sound. Can also be treated as a
% % function a(h,s), hence used in the above equation to find h*

% Define residual function
equation = @(h_star) log(2*(h0 - h_star)) - (2*sos_interp([log(s0), log(h_star)]));

% Define bounds for bisection
h_lower = 0;
h_upper = h0;

% Tolerance for the solution
tol = 1e-9;

% bisection method
while abs(h_upper - h_lower) > tol
    h_mid = (h_lower + h_upper) / 2;
    if equation(h_mid) * equation(h_lower) < 0
        h_upper = h_mid;
    else
        h_lower = h_mid;
    end
end

% Compute h_star as the midpoint of the final interval
h_star = (h_lower + h_upper) / 2;

% s_star remains constant
s_star = s0;

%final residual
residual = sqrt(2*(h0 - h_star)) - exp(sos_interp([log(s_star), log(h_star)]));
fprintf('\n\nFinal Residual for h*: %e\n', residual);
end