function [a_star, rho_star, F_rho_star_a_star] = compute_rho_star_astar_Fstar...
    (s0, h_star, rbf_interpolator_rho, rbf_interpolator_speed)

% this function finds the density and speed of sound at the throat

rho_star = exp(rbf_interpolator_rho(log([s0, h_star])));
a_star = exp(rbf_interpolator_speed(log([s0, h_star])));

F_rho_star_a_star = a_star*rho_star;
end