close all
clear all

STL_ac = 'fuselage_A320.stl';

triangles = read_binary_stl_file(STL_ac)
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates Geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% XFLR5_file = strcat('fuse_pepi_desplazado.txt');
% XFLR5_file = strcat('Fuselage_DATA_scaled.txt');
% STL_file = strcat('Fus_ProVant_v4.0.stl');

%--------------------------- CALCULATES GEOMETRY DATA ---------------------------------
% Data to read XFLR5 files
nSections = 9;
nPoints = 50;
lecture_Geo = 5;
STL_PLOT = 1;


% % Scale from fuselage in XFLR5 to CAD
% l_fus1 = 0.8799
% l_fus2 = Geo_tier.l_fus
% Ratio_SCALE = l_fus2/l_fus1
% pause
% Reads to extract length of fuselage
% assumes Ratio_SCALE = 1
Ratio_SCALE_l_fus_tmp = 1;
Ratio_SCALE_w_fus_tmp = 1;
Ratio_SCALE_h_fus_tmp = 1;
Ratio_SCALE_tmp.Ratio_SCALE_l_fus = Ratio_SCALE_l_fus_tmp;
Ratio_SCALE_tmp.Ratio_SCALE_w_fus = Ratio_SCALE_w_fus_tmp;
Ratio_SCALE_tmp.Ratio_SCALE_h_fus = Ratio_SCALE_h_fus_tmp;


Ratio_SCALE_l_fus = Ratio_SCALE_tmp.Ratio_SCALE_l_fus;
Ratio_SCALE_w_fus = Ratio_SCALE_tmp.Ratio_SCALE_w_fus;
Ratio_SCALE_h_fus = Ratio_SCALE_tmp.Ratio_SCALE_h_fus;

n_ele_geo = nPoints; % elementos que se han de leer
n_sections_body = nSections; % Elemtos del fuselaje

for i=1:n_sections_body
    lecture_Geo = lecture_Geo;
    % Extracción de datos VLM
    file_id=fopen('Fuselage_DATA_scaled.txt');
    aux = textscan(file_id,'%n%n%n','Headerlines',lecture_Geo);
    %      size(aux{i})
    %      pause
    Data_P(i).raw_data= aux;
    Data_P(i).label = 'Fuselage_DATA_scaled.txt';
    Data_P(i).x = aux{1};
    Data_P(i).y = aux{2};
    Data_P(i).z = aux{3};
    x(:,i)      = aux{1};
    y(:,i)      = aux{2};
    z(:,i)      = aux{3};
    x_max(i) = max(aux{1});
    y_max(i) = max(aux{2});
    z_max(i) = max(aux{3});
    fclose(file_id);
    lecture_Geo = lecture_Geo + n_ele_geo  + 2;
end

plotData{1}  = [x(:,end:(-1):1),x];
plotData{2}  = [-y(:,end:(-1):1), y];
plotData{3} = [z(:,end:(-1):1), z];

size(plotData{1})
size(plotData{2})
size(plotData{3})

figure(2)
subplot(1,3,1)
plot(plotData{1},plotData{2})
grid on
subplot(1,3,2)
plot(plotData{1},plotData{3})
grid on
subplot(1,3,3)
plot(plotData{2},plotData{3})
grid on

figure(3)
mesh(plotData{1},plotData{2},plotData{3});

% Loads STL file into a mesh
model = stl2matlab(STL_ac);
ScaleFactor = 1/1000; % CAD designs are in mm, and changing to m
plotData{1} = model{1}*ScaleFactor;
plotData{2} = model{2}*ScaleFactor;
plotData{3} = model{3}*ScaleFactor;

size(plotData{1})
size(plotData{2})
size(plotData{3})

figure(4)
subplot(1,3,1)
plot(plotData{1},plotData{2})
grid on
subplot(1,3,2)
plot(plotData{1},plotData{3})
grid on
subplot(1,3,3)
plot(plotData{2},plotData{3})
grid on

figure(5)
mesh(plotData{1},plotData{2},plotData{3});

% Rewrites data in the format of the XFLR5 lecture files in order to be used the
x_rec(:,1) = plotData{1}(:,3);
x_rec(:,2) = plotData{1}(:,4);
y_rec(:,1) = plotData{2}(:,3);
y_rec(:,2) = plotData{2}(:,4);
z_rec(:,1) = plotData{3}(:,3);
z_rec(:,2) = plotData{3}(:,4);

figure(6)
subplot(1,3,1)
plot(plotData{1},plotData{2})
grid on
subplot(1,3,2)
plot(plotData{1},plotData{3})
grid on
subplot(1,3,3)
plot(plotData{2},plotData{3})
grid on

figure(7)
subplot(1,3,1)
plot(plotData{1}(:,1000),plotData{2}(:,1000))
grid on
subplot(1,3,2)
plot(plotData{1}(:,1000),plotData{3}(:,1000))
grid on
subplot(1,3,3)
plot(plotData{2}(:,1000),plotData{3}(:,1000))
grid on

figure(8)
mesh(plotData{1}(:,1000),plotData{2}(:,1000),plotData{3}(:,1000));
            