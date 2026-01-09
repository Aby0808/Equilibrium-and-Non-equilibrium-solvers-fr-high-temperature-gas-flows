function [filtered_reactions, constants] = FUN_get_reactions(file, elements)
    % This function gets the reactions required from the reactions file 
    % based on the elements present in the reaction.

    % Open the Reactions.dat file and read lines
    fid = fopen(file, 'r');
    if fid == -1
        error('Cannot open Reactions.dat file');
    end

    % Initialize a cell array to hold the filtered reactions and constants
    filtered_reactions = {};
    constants = [];
    
    line = fgetl(fid);
    
    % Process each line until the end of the file
    while ischar(line)
        % Skip empty lines
        if isempty(line)
            line = fgetl(fid); 
            continue;
        end

        % If it's a reaction line (contains '<->'), process it
        if contains(line, '<->')
            % Split the reaction into reactants and products
            parts = strsplit(line, '<->');

            if length(parts) ~= 2
                line = fgetl(fid); % Skip malformed lines
                continue;
            end

            % Combine reactants and products to check for elements
            all_species = strjoin(parts, ' ');  % Join with space for better separation
            species = regexp(all_species, '[^\s]+', 'match');  % Extract all species

            % Check if all species in the reaction are in the specified elements set
            valid_reaction = true;
            for i = 1:length(species)
                % Check if each species is in the elements list
                if ~is_member_of_elements(species{i}, elements)
                    valid_reaction = false; % Found an element not in the set
                    break;
                end
            end

            if valid_reaction
                % Parse the constants in the next line
                line = fgetl(fid);
                constants_line = str2num(line); % Convert the constants to a numeric array
                if ~isempty(constants_line) % Ensure it parsed correctly
                    % Append the valid reaction and its constants
                    filtered_reactions{end + 1} = parts{1} + " <-> " + parts{2}; % Reaction string
                    constants = [constants; constants_line]; % Append constants
                end
            end
        end
        line = fgetl(fid); % Read the next line
    end

    fclose(fid); % Close the file
end

function isValid = is_member_of_elements(species, elements)
    % Check if species is a valid element or a valid multi-atom species
    if ismember(species, elements)
        isValid = true;
        return;
    end

    % Extract individual elements from the species string
    individual_elements = regexp(species, '[A-Z][a-z]*', 'match');
    
    % Check if all individual elements are in the elements list
    isValid = all(ismember(individual_elements, elements));
end
