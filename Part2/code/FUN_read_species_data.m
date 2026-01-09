function data = FUN_read_species_data(filename)
    % Initialize arrays and cell array for species names
    species_names = {}; % Cell array to store species names
    masses = [];        % Array to store masses
    bond_lengths = [];  % Array to store bond lengths
    frequencies = [];   % Array to store frequencies
    symmetries = [];    % Array to store symmetries

    % Open the file for reading
    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open file %s', filename);
    end

    % Read data
    while ~feof(fid)
        line = fgetl(fid);
        
        % Store the original species name
        species_names{end + 1} = line; %#ok<AGROW>
        
        % Sanitize species name for use as field name
        sanitized_line = matlab.lang.makeValidName(line);
        
        if startsWith(line, 'N') || startsWith(line, 'O') || startsWith(line, 'e')
            data_values = str2num(fgetl(fid));
            masses(1,end + 1) = data_values(1); % Store masses
            masses(2,end) = data_values(2); % Store masses
            bond_lengths(end + 1) = data_values(3); % Store bond length
            frequencies(end + 1) = data_values(4);  % Store frequency
            symmetries(end + 1) = data_values(5);   % Store symmetry
        else
            mass = str2double(fgetl(fid));
            masses(end + 1) = mass; % Store mass as a single value
        end
    end

    fclose(fid);
    
    % Create a single struct to hold all data
    data = struct();
    data.mass = masses;           % Cell array of masses
    data.bond_length = bond_lengths; % 1D array for bond lengths
    data.frequency = frequencies;   % 1D array for frequencies
    data.symmetry_factor = symmetries;     % 1D array for symmetries
    data.species_names = species_names; % Field for species names
end
