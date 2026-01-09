% Find Closest Index
function index = find_closest_index(array, value)
    % Find the index of the entry in 'array' that is closest to the 'value'
    % Ensure the array is in MATLAB format
    array = double(array); % Convert array to double if necessary
    
    % Find the index of the smallest difference
    [~, index] = min(abs(array - value));
end
