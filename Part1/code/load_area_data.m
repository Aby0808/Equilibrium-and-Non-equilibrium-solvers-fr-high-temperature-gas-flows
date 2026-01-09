function area = load_area_data(file_path)
    % Load the area variation data from file
    % Initialize an empty matrix to store the data
    fprintf('\nreading nozzle area data....')
    data = [];
    
    % Open the file
    fileID = fopen(file_path, 'r');
    if fileID == -1
        error('Could not open file: %s', file_path);
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
    expectedColumns = 2;
    if size(data, 2) ~= expectedColumns
        error('Unexpected data format: Expected %d columns, but got %d', expectedColumns, size(data, 2));
    end
    
    % Convert to table for easier usage
    area = array2table(data, 'VariableNames', {'x', 'A'});
    fprintf('\n..done!')
end
