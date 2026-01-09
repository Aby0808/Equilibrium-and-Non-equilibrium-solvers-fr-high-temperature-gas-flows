clear
close all

%author : Abhyudaya Singh
%NetID  : singh124@illinois.edu

%this program computes the properties inside a C-D nozzle with the flow in
%thermodynamic equilibrium as well as frozen flow

%% file read

% Load thermodynamic data from a file (contains enthalpy, entropy, density, etc.)
filename = 'data/output.dat';
[data, Enthalpy, Entropy, rho, speed_of_sound, Pressure, Temperature] = load_thermodynamic_data(filename);

% Load the area variation data: File containing the x and A(x) data
file_path = 'data/area.dat';
Area = load_area_data(file_path);

%% creating interpolators

% Create Radial Basis Function (RBF) interpolators for various thermodynamic properties
[rbf_rho, rbf_speed, rbf_pressure, rbf_temperature] = ...
    construct_rbf_interpolators(Enthalpy, Entropy, rho, speed_of_sound, Pressure, Temperature);

% Compute reservoir conditions based on the provided pressure and temperature
% These conditions define the thermodynamic state at the reservoir (upstream of the nozzle)
[enthalpy_interpolator, entropy_interpolator] = create_reservoir_interpolator(data);

%% reservoir conditions
% Example reservoir conditions: given pressure p and temperature T, get h and s
p0 = 5000000; % Reservoir pressure (Pa)
T0 = 4500;    % Reservoir temperature (K)
[h0, s0] = get_reservoir_h_and_s(p0, T0, enthalpy_interpolator, entropy_interpolator);

%% throat conditions

% Find the nozzle throat location (x*) and the corresponding minimum area (A*)
[x_star, A_star] = find_astar(Area);
index_star = find_closest_index(Area.A, A_star);

% Compute h* and s* at the throat using the speed of sound interpolator
[h_star, s_star] = compute_hstar_sstar(s0, h0, rbf_speed);

% Compute a*, rho* and F(h0, s0)
% [sound_star, rho_star, F_rho_a_star] = compute_rho_star_astar_Fstar(s_star, h_star, rbf_rho, rbf_speed);
[sound_star, rho_star, F_rho_a_star] = compute_rho_star_astar_Fstar(s_star, h_star, rbf_rho, rbf_speed);

%% solving the LTE flow in CD nozzle

A_x = Area.A;
result = process_nozzle_indirect_method(s0, h0, h_star, Area, A_x, A_star, F_rho_a_star, index_star, ...
    rbf_rho, rbf_speed, rbf_pressure, rbf_temperature);

% Extract results from the indirect method
enthalpy_values = result(:,1);
velocity_values = result(:,2);
density_values = result(:,3);
pressure_values = result(:,4);
temperature_values = result(:,5);
mach_values = result(:,6);
x_positions = result(:,7);

%plotting
figure(1)
plot(x_positions(116:end),enthalpy_values(116:end),'LineWidth',2)
hold on
grid on
xlabel('x (m)')
ylabel('Enthalpy (J/kg)')

figure(2)
plot(x_positions(116:end),velocity_values(116:end),'LineWidth',2)
hold on
grid on
xlabel('x (m)')
ylabel('Velocity (m/s)')

figure(3)
plot(x_positions(116:end),pressure_values(116:end),'LineWidth',2)
hold on
grid on
xlabel('x (m)')
ylabel('Pressure (Pa)')

figure(4)
plot(x_positions(116:end),temperature_values(116:end),'LineWidth',2)
hold on
grid on
xlabel('x (m)')
ylabel('Temperature (K)')

figure(5)
plot(x_positions(116:end),density_values(116:end),'LineWidth',2)
hold on
grid on
xlabel('x (m)')
ylabel('Density (Kg/m3)')

figure(6)
plot(x_positions(116:end),mach_values(116:end),'LineWidth',2)
hold on
grid on
xlabel('x (m)')
ylabel('Mach number')


%% Frozen flow

% Constants for frozen flow
R = 8.3144598 / 0.02672963120279829;
gamma = 1.242430995249157;

% solving for frozen flow
result_fr = process_nozzle_perfect_gas(gamma, R, p0, T0, Area, A_x, A_star, index_star);

enthalpy_values_fr = result_fr.enthalpy;
velocity_values_fr = result_fr.velocity;
density_values_fr = result_fr.density;
pressure_values_fr = result_fr.pressure;
temperature_values_fr = result_fr.temperature;
mach_values_fr = result_fr.mach;
x_positions_fr = result_fr.x;

%plotting
figure(1)
plot(x_positions_fr(1:end-1),enthalpy_values_fr(1:end-1),'LineWidth',2)
plot([x_star,x_star],[0,7.5*10^6],'Color','k')
legend('LTE flow','Frozen flow','throat')

figure(2)
plot(x_positions_fr(1:end-1),velocity_values_fr(1:end-1),'LineWidth',2)
plot([x_star,x_star],[0,3000],'Color','k')
legend('LTE flow','Frozen flow','throat')

figure(3)
plot(x_positions_fr(1:end-1),pressure_values_fr(1:end-1),'LineWidth',2)
plot([x_star,x_star],[0,5*10^6],'Color','k')
legend('LTE flow','Frozen flow','throat')

figure(4)
plot(x_positions_fr(1:end-1),temperature_values_fr(1:end-1),'LineWidth',2)
plot([x_star,x_star],[0,5000],'Color','k')
legend('LTE flow','Frozen flow','throat')

figure(5)
plot(x_positions_fr(1:end-1),density_values_fr(1:end-1),'LineWidth',2)
plot([x_star,x_star],[0,4],'Color','k')
legend('LTE flow','Frozen flow','throat')

figure(6)
plot(x_positions_fr(1:end-1),mach_values_fr(1:end-1),'LineWidth',2)
plot([x_star,x_star],[0,3],'Color','k')
legend('LTE flow','Frozen flow','throat')