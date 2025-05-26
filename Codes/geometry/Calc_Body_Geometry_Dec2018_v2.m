function [Body_Geo, plotData] = Calc_Body_Geometry_Dec2018_v2(XFLR5_file,STL_file,Geo_tier,...
    nPoints, nSections,conversion_UNITS,lecture_Geo,plotData,STLnSections,STLnPoints)

R2D = conversion_UNITS.R2D ;

% lecture_Geo = 5; % lineas a partir desde donde empieza a leer
n_ele_geo = nPoints; % elementos que se han de leer
n_sections_body = nSections; % Elemtos del fuselaje

Sref = Geo_tier.S_w1;

for i=1:n_sections_body
    lecture_Geo = lecture_Geo;
    % Extracción de datos VLM
    file_id=fopen(XFLR5_file);
    aux = textscan(file_id,'%n%n%n','Headerlines',lecture_Geo);
    Data_P(i).raw_data= aux;
    Data_P(i).label = XFLR5_file;
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

% % Loads STL file into a mesh
% model = stl2matlab('ProVant.stl');
% ScaleFactor = 1/1000; % CAD designs are in mm, and changing to m
% plotData{1} = model{1}*ScaleFactor;
% plotData{2} = model{2}*ScaleFactor;
% plotData{3} = model{3}*ScaleFactor;
% size(plotData{1})

% plotData{1}  = [x(:,end:(-1):1),x];
% plotData{2}  = [-y(:,end:(-1):1), y];
% plotData{3} = [z(:,end:(-1):1), z];

% Rewrites data in the format of the XFLR5 lecture files in order to be used in the codes 
x_rec(:,1) = plotData{1}(:,3);
x_rec(:,2) = plotData{1}(:,4);
y_rec(:,1) = plotData{2}(:,3);
y_rec(:,2) = plotData{2}(:,4);
z_rec(:,1) = plotData{3}(:,3);
z_rec(:,2) = plotData{3}(:,4);

% Calculo de los perímetros de cada sección
for i=1:(n_sections_body)
    Per(i).y = [Data_P(i).y; -Data_P(i).y;Data_P(i).y(1)];
    Per(i).z = [Data_P(i).z; fliplr(Data_P(i).z')';Data_P(i).z(1)];
end

% Calculo de los perímetros de cada sección
for i=1:(STLnSections)
    Per(i).y = [y_rec(i,:)'; -y_rec(i,:)'; y_rec(1)];
    Per(i).z = [z_rec(i,:)'; fliplr(z_rec(i,:))';z_rec(1)];
    x_max(i) = max(x_rec(i,:)');
    y_max(i) = max(y_rec(i,:)');
    z_max(i) = max(z_rec(i,:)');
end

x = x_rec;
y = y_rec;
z = z_rec;

n_sections_body = STLnSections;
% plotData{1}  = [x(:,end:(-1):1),x];
% plotData{2}  = [-y(:,end:(-1):1), y];
% plotData{3} = [z(:,end:(-1):1), z];

% for i=1:(n_sections_body-2)
for i=2:(n_sections_body)   
%     [geom,iner,cpmo] = polygeom(Per(i+1).y,Per(i+1).z);
    [geom,iner,cpmo] = polygeom(Per(i).y,Per(i).z);
    Area(i) = geom(1);
    X_cen(i) = geom(2);
    Y_cen(i) = geom(3);
    Z_cen(i) = geom(3); % assumes y and Z are the same
    Perimeter(i) = geom(4);
end

Body_Geo.S_x    = Area;
Body_Geo.x      = x_max;
% Body_Geo.l      = max(Body_Geo.x);

Delta_h             = Body_Geo.x(2:end) - Body_Geo.x(1:(end-1));
Delta_Zcent         = Z_cen(2:end)-Z_cen(1:(end-1));
size(Delta_h)
size(Delta_Zcent)
alpha               = -atan(Delta_h./Delta_Zcent);
Body_Geo.alpha_x    = [alpha,alpha(end)];

% Calculo del volumen (aproximación)
l=1;
for i=3:n_sections_body-1
    Delta_h = x_max(i) - x_max(i-1);
    Surface(l) = Delta_h*(Perimeter(i-1) + Perimeter(i-2))/2;
    l = l+1;
end

% Surface total geometría
Surf_TOT = sum(Surface);
Body_Geo.Surf_TOT = Surf_TOT;

% Calculo de las áreas de cada sección
for i=1:n_sections_body
    Area_body(i) = 2*polyarea(Data_P(i).y,Data_P(i).z);
%     dSdX(i) = diff(Area_body(i))/diff(x_max(i));
end
    
% Body_Geo.Area_body = Area_body;
Body_Geo.Area_body = Area_body;

% X-positions of fuselage
Body_Geo.x_Area_body = x_max;

% Calculo de los perímetros de cada sección
for i=1:(n_sections_body)
    Height(i) = max(Per(i).z) - min(Per(i).z);
    Bisectriz_z(i) = Height(i)/2;
    Centroid_z(i) =  max(Per(i).z) - Bisectriz_z(i);
end

Angle_fus = gradient(Centroid_z,x_max)*R2D;
x_interp = linspace(x_max(1),x_max(end),200);
Centroid_z_interp = interp1(x_max,Centroid_z,x_interp,'spline');
Angle_fus_interp = gradient(Centroid_z_interp,x_interp)*R2D;

Body_Geo.Angle_fus = Angle_fus;
Body_Geo.Angle_fus_interp = Angle_fus_interp;
Body_Geo.x_interp = x_interp;

% Calculo del ángulo de cada sección
for i=2:n_sections_body
    Delta_h(i-1) = x_max(i) - x_max(i-1);
%     Angle_fus(i-1) = atan(Centroid_z(i)/Delta_h(i-1));
end

Body_Geo.Height = Height;
Body_Geo.Bisectriz_z = Bisectriz_z;
Body_Geo.Centroid_z = Centroid_z;
Body_Geo.Angle_fus = Angle_fus;

% Calculo del volumen (aproximación)
l=1;
for i=2:n_sections_body
    Delta_h = x_max(i) - x_max(i-1);
    Vol(l) = Delta_h*(Area_body(i) + Area_body(i-1))/2;
    l = l+1;
end

% Volumen total geometría
Vol_TOT = sum(Vol);
Body_Geo.Vol_TOT = Vol_TOT;

% Max area and location
[Area_b_max i_max] = max(Area_body);
x_Area_b_max = x_max(i_max);
w_Area_b_max = 2*y_max(i_max);
h_Area_b_max = 2*z_max(i_max);

% Calculo de altitud en sección 1/4 y 3/4 fuselaje
x_Area_b_1_4 = x_max(end)/4;
x_Area_b_3_4 = x_max(end)*(3/4);

% Calculo de altitud en sección 1/4 y 3/4 fuselaje
z_Area_b_1_4 = 2*interp1(x_max,z_max,x_Area_b_1_4,'spline');
z_Area_b_3_4 = 2*interp1(x_max,z_max,x_Area_b_3_4,'spline');

Body_Geo.z_Area_b_1_4 = z_Area_b_1_4;
Body_Geo.z_Area_b_3_4 = z_Area_b_3_4;

Body_Geo.x_Area_b_1_4 = x_Area_b_1_4;
Body_Geo.x_Area_b_3_4 = x_Area_b_3_4;

Body_Geo.Area_b_max = Area_b_max;
Body_Geo.x_Area_b_max = x_Area_b_max;
Body_Geo.w_Area_b_max = w_Area_b_max;
Body_Geo.h_Area_b_max = h_Area_b_max;

% Calculo del Area proyectada en plano x-y
Area_top  = 2*polyarea(x_max,y_max);
% Calculo del Area proyectada en plano x-z
Area_side = 2*polyarea(x_max,z_max);

Body_Geo.Area_top = Area_top;
Body_Geo.Area_side = Area_side;

Body_Geo.l_fus = x_max(end);

length_x_position = x_max;
width_x_position = 2*y_max;
height_x_position = 2*z_max;

Body_Geo.length_x_position = length_x_position;
Body_Geo.width_x_position = width_x_position;
Body_Geo.height_x_position = height_x_position;

% OLD CODES
% Calculo del volumen (aproximación)
Vol = Delta_h.*(Body_Geo.S_x(2:end) + Body_Geo.S_x(1:(end-1)))/2;

% Volumen total geometría
dSdX                = (Body_Geo.S_x(2) - Body_Geo.S_x(1))/(Body_Geo.x(2) - Body_Geo.x(1));
dSdX                = [dSdX,(Body_Geo.S_x(3:end) - Body_Geo.S_x(1:(end-2)))./(Body_Geo.x(3:end) - Body_Geo.x(1:(end-2)))];
dSdX                = [dSdX, (Body_Geo.S_x(end) - Body_Geo.S_x(end-1))/(Body_Geo.x(end) - Body_Geo.x(end-1))];
Body_Geo.dSdX       = dSdX;

i_SfrontMax = 1;
while dSdX(i_SfrontMax) > 0
    i_SfrontMax=i_SfrontMax+1;
end

Body_Geo.S_front = Body_Geo.S_x(i_SfrontMax);
Body_Geo.x_S_max  = x_max(i_SfrontMax);
w_Area_b_max    = 2*y_max(i_SfrontMax);
h_Area_b_max    = 2*z_max(i_SfrontMax);
Body_Geo.w_Area_b_max = w_Area_b_max;
Body_Geo.h_Area_b_max = h_Area_b_max;

% Calculo del Area proyectada en plano x-y
Body_Geo.S_top  = 2*polyarea(x_max,y_max);
% Calculo del Area proyectada en plano x-z (se asume seccion simetrica)
Body_Geo.S_side = 2*polyarea(x_max,z_max);

%length_x_position = x_max;
% Width
Body_Geo.W_x = 2*y_max;
Body_Geo.w_fus   = max(Body_Geo.W_x);
% Height
Body_Geo.H_x = 2*x_max; 
Body_Geo.h_fus   = max(Body_Geo.H_x);
% Diameter
Body_Geo.D_x = (Body_Geo.H_x + Body_Geo.W_x)/2;
Body_Geo.d_fus   = max(Body_Geo.D_x);

Body_Geo.x_max   = x_max;
Body_Geo.y_max   = y_max;
Body_Geo.z_max   = z_max;

%% PROPIEDADES AERODINÁMICAS
% Cálculo de la pendiente de sustentación del fuselaje
f_k2_k1         = k2_k1_calc(Body_Geo.l_fus, Geo_tier.d_fus);      % Fuselage apparent mass coefficient. Pamadi, Figure 3.6
Body_Geo.CLa_fus    = 2*f_k2_k1*(Body_Geo.S_front/Sref);         % Pendiente de sustentación del morro aislado
end
