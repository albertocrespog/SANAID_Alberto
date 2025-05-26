function [slice_height, z_range] = calculate_slice_height(triangles, num_slices)
    % Calculate the range in the Z direction
    z_min = min(triangles(:, [3, 6, 9]), [], 'all');
    z_max = max(triangles(:, [3, 6, 9]), [], 'all');
    z_range = z_max - z_min;

    % Calculate the slice height to achieve the desired number of slices
    slice_height = z_range / num_slices;
end