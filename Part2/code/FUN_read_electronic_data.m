function data = FUN_read_electronic_data(filename)
    % Initialize the struct
    data.g_elec = [];
    data.energies = [];
    data.species_names = {}; % Cell array to store species names

    % Open the file for reading
    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open file %s', filename);
    end

    % Read data
    while ~feof(fid)
        line = fgetl(fid); % Read the current line

        if ~ischar(line) || isempty(line)
            continue; % Skip empty lines or if line is not a char
        end
        
        % Store the current species name
        current_species = line;
        data.species_names{end + 1} = current_species; % Append the species name
        
        % Initialize arrays for the current species
        current_g_elec = [];
        current_energies = [];

        while true
            line = fgetl(fid); % Read the next line
            
            if ~ischar(line) || isempty(line)
                break; % Stop if we hit an empty line or end of file
            end
            
            % Check if the line starts with a letter indicating a new species
            if startsWith(line, 'N') || startsWith(line, 'O') || startsWith(line, 'e')
                % Move the file pointer back for the next iteration
                fseek(fid, -length(line) - 2, 'cof'); 
                break; % Stop reading if we hit a new species
            end
            
            % Read the data values
            data_values = sscanf(line, '%f'); % Use sscanf to read numbers
            
            % Ensure we have two values
            if numel(data_values) == 2
                current_g_elec(end + 1) = data_values(1); % First value is g
                current_energies(end + 1) = data_values(2); % Second value is energy in cm^-1
            end
        end
        
        % Concatenate current values to main arrays
        data.g_elec = [data.g_elec; current_g_elec]; % Append current g_elec
        data.energies = [data.energies; current_energies]; % Append current energies
    end

    fclose(fid);
end
