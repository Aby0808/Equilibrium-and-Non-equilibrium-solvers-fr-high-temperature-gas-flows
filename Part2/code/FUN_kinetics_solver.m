function [om_s] = FUN_kinetics_solver(Tt,Tv, rho, y)

% this function computes the forward and backward reaction rates and
% calculates the mass rate of production of speceis

global ARR_CONST SP_AVAIL SP_MAT IND_MAT ST_COEFF_MAT

T = sqrt(Tt*Tv);
%% initializing necessary variables
num_species = length(SP_AVAIL); % num of species available [N2,N]
num_react = length(SP_MAT);  % Number of reactions
Kfb = zeros(num_react, 2);      % Preallocate K array
Kc = FUN_get_Kp_eq_stat(Tt);     % equilibrium constant from stat thermo

Rfr = ones(num_react,1);
Rbr = ones(num_react,1);
om_s = zeros(num_species,1);

%% calculating reaction rates

% Loop through each reaction
for i = 1:num_react
    nb_sp_react = length(SP_MAT{i});  % Number of species in the current reaction
    Kfb(i,1) = (ARR_CONST(i,1)* 10^ARR_CONST(i,2))*(T^ARR_CONST(i,3))*exp(-ARR_CONST(i,4)/T); %forward rate
    Kfb(i,2) = Kfb(i,1)/Kc;                                                               % backward rate
    %Kfb(i,2) = (ARR_CONST(i,5)^ARR_CONST(i,6))*(T^ARR_CONST(i,7));
    % Loop through each species in the current reaction
    for j = 1:nb_sp_react
        current_species = SP_MAT{i}(j);
        Ms = FUN_get_thermo_prop(current_species, T, 'mol_wt'); %converting from kg/mol to kg/kg-mol
        if strcmp(current_species,'N2')
            rhos = rho*y(1);
        elseif strcmp(current_species,'N')
            rhos = rho*y(2);
        else
            error('invalid species %s in reaction',current_species)
        end
        if strcmp(IND_MAT(i,j),'r')
            Rfr(i) = Rfr(i)*0.001*(rhos/Ms)^ST_COEFF_MAT(i,j);
        elseif strcmp(IND_MAT(i,j),'p')
            Rbr(i) = Rbr(i)*0.001*(rhos/Ms)^ST_COEFF_MAT(i,j);
        end
    end
    % Rfr(i) = 0.001*Kfb(i,1)*Rfr(i);
    % Rbr(i) = 0.001*Kfb(i,2)*Rbr(i);
    Rfr(i) = Kfb(i,1)*Rfr(i);
    Rbr(i) = Kfb(i,2)*Rbr(i);
end

%% calculating the mass rate of production of species

% looping through each species available
for k = 1:num_species
    sum = 0;
    for i = 1:num_react
        nb_sp_react = length(SP_MAT{i});  % Number of species in the current reaction
        % Kfb(i,1) = ARR_CONST(i,1)*(T^ARR_CONST(i,2))*exp(-ARR_CONST(i,3)/T); %calculating rate constants
        % Kfb(i,2) = Kc/Kfb(i,1);
        % Loop through each species in the current reaction
        alpha = 0; beta = 0;
        for j = 1:nb_sp_react
            current_species = SP_MAT{i}(j);

            if strcmp(current_species,SP_AVAIL(k))
                if strcmp(IND_MAT(i,j),'r')
                    alpha = alpha + ST_COEFF_MAT(i,j);
                elseif strcmp(IND_MAT(i,j),'p')
                    beta = beta + ST_COEFF_MAT(i,j);
                end
            end
        end
        sum = sum + (beta-alpha)*(Rfr(i)-Rbr(i));
    end
    Ms = FUN_get_thermo_prop(SP_AVAIL(k), T, 'mol_wt'); %converting from kg/mol to kg/kg-mol
    om_s(k) = Ms*sum;
end
end
