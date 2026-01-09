
% Find A* and its location
function [x_star, A_star] = find_astar(area)
    % Find the location and value of A*
    [A_star, idx_min] = min(area.A); % Find the minimum value and its index
    x_star = area.x(idx_min);        % Get the corresponding x value
end
