function [points, triangles, tri_norms] = import_stl_ascii(filename, mode)
% Import ASCII STL files
fid = fopen(filename, 'r');
fmt = '%*s %*s %f %f %f \r\n %*s %*s \r\n %*s %f %f %f \r\n %*s %f %f %f \r\n %*s %f %f %f \r\n %*s \r\n %*s \r\n';
C = textscan(fid, fmt, 'HeaderLines', 1);
fclose(fid);

% Extract normal vectors and vertices
tri_norms = cell2mat(C(1:3));
tri_norms = tri_norms(1:end-1, :); % Remove junk from last line

v1 = cell2mat(C(4:6));
v2 = cell2mat(C(7:9));
v3 = cell2mat(C(10:12));

if isnan(C{4}(end))
    v1 = v1(1:end-1, :); % Remove junk from last line
    v2 = v2(1:end-1, :);
    v3 = v3(1:end-1, :);
end

v_temp = [v1 v2 v3]';
v = zeros(3, numel(v_temp) / 3);
v(:) = v_temp(:);
v = v';

switch mode
    case 1
        [points, triangles] = fv2pt(v, length(v) / 3); % Get points and triangles
        tri_norms = tri_norms;
    case 2
        points = v;
        triangles = tri_norms;
end
end 