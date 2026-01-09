function [data, Enthalpy, Entropy, rho, speed_of_sound, Pressure, Temperature] = load_thermodynamic_data(filename)
    % Load thermodynamic data from the specified file
    % Initialize an empty matrix to store the data
    fprintf('\nreading thermo data...')
    data = [];
    
    % Open the file
    fileID = fopen(filename, 'r');
    if fileID == -1
        error('Could not open file: %s', filename);
    end
    
    % Read the file line by line
    line = fgetl(fileID);
    while ischar(line)
        % Ignore comment lines starting with '#'
        if startsWith(strtrim(line), '#')
            line = fgetl(fileID);
            continue;
        end
        
        % Parse numerical values from the line
        lineData = sscanf(line, '%f');
        if ~isempty(lineData)
            data = [data; lineData']; % Append row-wise
        end
        
        % Read the next line
        line = fgetl(fileID);
    end
    
    % Close the file
    fclose(fileID);
    
    % Check if the data has the expected number of columns
    expectedColumns = 7;
    if size(data, 2) ~= expectedColumns
        error('Unexpected data format: Expected %d columns, but got %d', expectedColumns, size(data, 2));
    end
    
    % Assign column names and extract individual columns
    Temperature    = data(:, 1); % Column 1
    Pressure       = data(:, 2); % Column 2
    rho            = data(:, 3); % Column 3
    MolarMass      = data(:, 4); % Column 4 (not used but available)
    Enthalpy       = data(:, 5); % Column 5
    Entropy        = data(:, 6); % Column 6
    speed_of_sound = data(:, 7); % Column 7
    fprintf('\n..done!')
end
