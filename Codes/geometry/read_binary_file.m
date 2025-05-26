function tr=read_binary_file(fn)
f = fopen(fn,'r');
rd = fread(f,inf,'uint8=>uint8');
numTriangles = typecast(rd(81:84),'uint32');
tr = zeros(numTriangles,12);
sh = reshape(rd(85:end),50,numTriangles);
tt = reshape(typecast(reshape(sh(1:48,1:numTriangles),1,48*numTriangles),'single'),12,numTriangles)';
tr(:,1:9) = tt(:,4:12);
tr(:,10:12) = tt(:,1:3);
end