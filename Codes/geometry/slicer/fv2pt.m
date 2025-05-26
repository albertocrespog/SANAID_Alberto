function [p, t] = fv2pt(v, fnum)
    % Get points and triangle indices given vertices and facet number
    c = size(v, 1);

    % Check if the number of vertices is divisible by 3
    if mod(c, 3) ~= 0
        warning('Number of vertices is not divisible by 3. Removing incomplete triangles.');
        c = floor(c / 3) * 3; % Adjust to the nearest multiple of 3
        v = v(1:c, :); % Truncate the extra vertices
    end

    % Recalculate fnum if inconsistent
    if c ~= 3 * fnum
        warning('Adjusting fnum to match the number of vertices.');
        fnum = c / 3;
    end

    % Triangles with vertex ID data
    t = zeros(3, fnum);
    t(:) = 1:c;

    % Keep unique points from vertices
    [p, ~, j] = unique(v, 'rows');
    t(:) = j(t(:));
    t = t';
end