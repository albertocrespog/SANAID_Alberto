function tr=read_ascii_file(fn)
fid = fopen(fn,'r');
tr=[];
tr_n=[];
tline='start';
while ~contains(tline,'endsolid')
    tline = fgetl(fid);
    if contains(tline, 'facet normal')
        tr_n = [tr_n; sscanf(tline, '%*s %*s %f %f %f')'];
    elseif contains(tline, 'vertex')
        tr = [tr; sscanf(tline, '%*s %f %f %f')'];
    end
end
tr=reshape(tr',9,[])';
tr=[tr,tr_n];
end