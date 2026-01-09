function [x,prim_var] = FUN_solve_1d_noneq_flow(x_span, y_initial, u_initial, T_initial, Tv_initial)
    % Solve 1D flow equations using ode15s

    prim_var0 = [y_initial, u_initial, T_initial, Tv_initial]; % Combine all initial values

    % Solve the ODE system using ode15s
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8); % Set tolerances for accuracy
    [x, prim_var] = ode15s(@ode_system, x_span, prim_var0, options); % solving the equations
end
