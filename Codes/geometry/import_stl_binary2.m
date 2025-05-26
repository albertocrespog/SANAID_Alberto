function [points, triangles, tri_norms] = import_stl_binary2(filename, mode)
% Import binary STL files
fid = fopen(filename, 'rb');

% Read header and number of facets
fseek(fid, 80, 'bof');
fnum = fread(fid, 1, 'uint32');

triangles = zeros(fnum, 9);
tri_norms = zeros(fnum, 3);

for i = 1:fnum
    tri_norms(i, :) = fread(fid, 3, 'float32'); % Normal vector
    vertex1 = fread(fid, 3, 'float32');
    vertex2 = fread(fid, 3, 'float32');
    vertex3 = fread(fid, 3, 'float32');
    fread(fid, 1, 'uint16'); % Skip attribute byte count
    
    triangles(i, :) = [vertex1', vertex2', vertex3'];
end

fclose(fid);

points = unique(reshape(triangles, [], 3), 'rows');
switch mode
    case 1
        [points, triangles] = fv2pt(points, fnum);
    case 2
        triangles = tri_norms;
end
end