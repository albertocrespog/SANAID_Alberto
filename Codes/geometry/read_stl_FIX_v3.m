function triangles = read_stl_FIX_v3(filename, scale_factor)
    % Detect the STL format and read accordingly
    % Supports both ASCII and binary STL files

    if nargin < 2
        scale_factor = 1; % Default scale factor
    end

    fid = fopen(filename, 'r');
    header = fread(fid, 80, 'uint8=>char')'; % Read the first 80 characters of the file
    fclose(fid);

    % Check if file is ASCII or binary
    if contains(header, 'solid') && ~contains(header, char(0))
        triangles = read_ascii_stl(filename, 1); % Call existing ASCII reader
    else
        triangles = read_binary_stl(filename); % Use a binary STL reader
    end

    % Apply scaling factor to the triangles
    triangles(:, 1:9) = triangles(:, 1:9) * scale_factor;
end



