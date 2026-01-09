clc
clear
close all

%author : Abhyudaya Singh
%NetID  : singh124@illinois.edu

%This program solves the non equilibrium flow behind a normal shock wave

%% setting global variables

% initializing global variables
fprintf('\ninitializing constants\n')
global MDOT C_LIGHT RU K_B H_PLANCK NA_AVO SPECIES_DATA SPEC_DATA SP_AVAIL ELEM_AVAIL REACTION ARR_CONST SP_MAT IND_MAT ST_COEFF_MAT
C_LIGHT = 2.998e8;
K_B = 1.380649e-23; % Boltzmann constant (J/K)
H_PLANCK = 6.62607015e-34; % Planck constant (J s)
NA_AVO = 6.023e23;      % Avogardo's number
RU = 8.314; % universal gas constant

SP_AVAIL = ["N2","N"];
ELEM_AVAIL = 'N';

%% file read
fprintf('\nReading necessary data\n')

% names of the relevant files
Species_data_file = 'data/Species_data.dat';
Spectroscopic_data_file = 'data/Atomic_spectra_database.dat';
Reactions_data_file = 'data/Reactions.dat';

% reading necessary data files
SPECIES_DATA=FUN_read_species_data(Species_data_file);          %species data required for calculating thermo properties through stat mech approach
SPEC_DATA=FUN_read_electronic_data(Spectroscopic_data_file);    %spectroscopic data for electronic partition functions 
[REACTION,ARR_CONST] = FUN_get_reactions(Reactions_data_file,ELEM_AVAIL); %reading the necessary reactions and constant for arrhenius eqn
[SP_MAT,IND_MAT,ST_COEFF_MAT] = FUN_get_Stoichiometric_Coeff();

fprintf('\n..done!\n\n\n')
%% Initialization 

fprintf('---- Initializing freestream-----\n')

% initializing the pre shock conditions

P1 = 5;     % Pre shock pressure
T1 = 200;   % Pre shock temperature
Tv1 = T1;   % Pre shock vibrational temperature
Vs = 7000;  % Shock speed 
x_max = 0.1;% Shock stand-off distance

% initializing the gas composition
[y,species] = initialize_composition();    % mass fractions

fprintf('\nShock-Tube simulation\n')
fprintf('=========================\n')
fprintf('Mixture: N2/N\n\n')

%% Compute initial conditions (Post shock)

fprintf(' Computing Rankine - Hugoniot jump relations for hot gas\n')
[P2, rho2, T2, V2] = compute_rh_jump(y,species,P1,T1,Vs);
fprintf('\nPost shock conditions (Hot gas): ')
fprintf('\nPressure: %d',P2)
fprintf('\nDensity: %d',rho2)
fprintf('\nTemperature: %d',T2)
fprintf('\nVelocity: %d\n',V2)

MDOT = V2*rho2;

%% Simulate the flowfield downstream of the shock
% constructing the partial differential equations and solving them to get downstream flowfield

xspan = [0,0.1]; % length of the domain
fprintf('\nsolving the post shock flow!\n')
[x,prim_var] = FUN_solve_1d_noneq_flow(xspan,y,V2,T2,Tv1);

%% extract data

% Extract results
num_species = length(SP_AVAIL);
y_post = prim_var(:, 1:num_species);   % Species mass fractions
u_post = prim_var(:, num_species + 1); % Velocity
T_post = prim_var(:, num_species + 2); % Translational temperature
Tv_post = prim_var(:, num_species + 3); % Vibrational temperature
rho_post = zeros(size(u_post));
P_post = zeros(size(u_post));
for i = 1:size(y_post,1)
    rho_post(i) = MDOT/u_post(i);
    R_sp = FUN_get_mix_thermo_prop(species,y_post(i,:),T_post(i),'R_sp');
    P_post(i) = rho_post(i)*R_sp*T_post(i);
end

fprintf('\ndone!\n')
%% plotting the results

% Plot the results
figure
plot(x, y_post);
grid on
title('Mass Fractions');
xlabel('x (m)');
ylabel('y_i');
legend(SP_AVAIL);

figure
plot(x, u_post);
grid on
title('Velocity');
xlabel('x (m)');
ylabel('u (m/s)');

figure
plot(x, T_post);
grid on
title('Translational temperature');
xlabel('x (m)');
ylabel('T (K)');

figure
semilogx(x, T_post);
hold on
semilogx(x, Tv_post);
grid on
title('Translational and Vibrational temperature');
xlabel('x (m)');
ylabel('T (K)');
legend('T','Tv')

figure
plot(x, rho_post);
grid on
title('Density');
xlabel('x (m)');
ylabel('rho (kg/m3)');

figure
plot(x, P_post);
grid on
title('Pressure');
xlabel('x (m)');
ylabel('P (Pa)');

%% plotting distribution function
FUN_plot_distributions('N2',x,T_post,Tv_post,[0,0.01,0.05,0.1]);
FUN_plot_distributions('N',x,T_post,Tv_post,[0,0.01,0.05,0.1]);