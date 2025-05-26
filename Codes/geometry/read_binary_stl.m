function triangles = read_binary_stl(filename)
    % Reads binary STL file
    fid = fopen(filename, 'rb');
    fseek(fid, 80, 'bof'); % Skip header
    n_triangles = fread(fid, 1, 'uint32');
    
    triangles = zeros(n_triangles, 12); % Initialize triangle matrix
    for i = 1:n_triangles
        normal = fread(fid, 3, 'float32'); % Read normal vector
        vertex1 = fread(fid, 3, 'float32'); % Read first vertex
        vertex2 = fread(fid, 3, 'float32'); % Read second vertex
        vertex3 = fread(fid, 3, 'float32'); % Read third vertex
        fread(fid, 1, 'uint16'); % Skip attribute byte count
        
        % Store data in the triangle matrix
        triangles(i, :) = [vertex1', vertex2', vertex3', normal'];
    end
    fclose(fid);
end