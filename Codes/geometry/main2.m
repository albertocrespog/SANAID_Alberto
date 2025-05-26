clc;clear;
filename='fuselage_Cessna208_v3_after.stl';
triangles = read_stl_file(filename);

tri=reshape(triangles(:,1:9)',[3,3*size(triangles,1)])';
T=reshape(1:length(tri),3,length(tri)/3)';
trimesh(T,tri(:,1),tri(:,2),tri(:,3),'FaceColor',[1 1 1],'EdgeColor','b');
axis equal off
hold on;