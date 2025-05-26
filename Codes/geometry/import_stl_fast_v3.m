function varargout=import_stl_fast_v3(filename,mode)
% Function to import STL files into MATLAB. Supports both ASCII and binary formats.
%--------------------------------------------------------------
% Inputs:       filename:   string 
%               mode:       integer 1 or 2, see outputs:
%
% Outputs:      mode 1: [points,triangles,tri_norms, stats]
%               mode 2: [vertices, tri norms, stats]
%
% Author: Eric Trautmann, updated to support binary STL files.
%--------------------------------------------------------------
% Ex: 
% filename = 'file.stl'     
% [p,t,tnorm,stats] = import_stl_fast_v2(filename,1)

if nargin<2
    mode=1; % Default value
end

if ~(mode==1 || mode==2)
    error('Invalid mode')
end

if nargout<3 && mode==1
    error('Invalid input number / mode setting')
end
if nargout>2 && mode==2
    error('Invalid input number / mode setting')
end

% Open file
fid = fopen(filename, 'r'); % Open the file
if fid == -1
    error('File could not be opened, check name or path.')
end

% Detect file format (ASCII or binary)
header = fread(fid, 80, 'uint8=>char')';
fseek(fid, 0, 'bof'); % Reset file position

if contains(header, 'solid') && ~contains(header, char(0))
    % ASCII STL
    fclose(fid);
    [varargout{1:nargout-1}] = import_stl_ascii(filename, mode);
else
    % Binary STL
    fclose(fid);
    [varargout{1:nargout-1}] = import_stl_binary2(filename, mode);
end

% Calculate stats (min, max, and deltas for each axis)
points = varargout{1};
min_vals = min(points);
max_vals = max(points);
deltas = points - min(points);

stats = struct('min', min_vals, 'max', max_vals, 'deltas', deltas);
varargout{end} = stats;
end





