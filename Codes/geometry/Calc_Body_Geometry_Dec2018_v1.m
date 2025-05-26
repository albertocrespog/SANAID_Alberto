function [Fig,Body_Geo,plotData] = Calc_Body_Geometry_Dec2018_v1(XFLR5_file,STL_files,Geo_tier,...
    nPoints, nSections,conversion_UNITS,lecture_Geo,Ratio_SCALE,SF,OUTPUT_read_XLSX,Fig)

R2D = conversion_UNITS.R2D;

Ratio_SCALE_l_fus = Ratio_SCALE.Ratio_SCALE_l_fus;
Ratio_SCALE_w_fus = Ratio_SCALE.Ratio_SCALE_w_fus;
Ratio_SCALE_h_fus = Ratio_SCALE.Ratio_SCALE_h_fus;

% lecture_Geo = 5; % lineas a partir desde donde empieza a leer
n_ele_geo = nPoints; % elementos que se han de leer
n_sections_body = nSections; % Elemtos del fuselaje

Sref = Geo_tier.S_w1;

%% Reads the STL file
% model = stl2matlab(OUTPUT_read_XLSX.Fuselage_flags.STL_fus)
SF_CAD = OUTPUT_read_XLSX.Fuselage_flags.SF_CAD;
l_fus = OUTPUT_read_XLSX.InputGeometry_Data_flags.l_fus;
% N_slices = OUTPUT_read_XLSX.InputGeometry_Data_flags.N_slices;
SF_CAD = OUTPUT_read_XLSX.Fuselage_flags.SF_CAD;
N_slices = 100;
% slice_height = (l_fus/N_slices)/SF_CAD;
slice_height = (l_fus/N_slices)/SF_CAD;
STL_fus = OUTPUT_read_XLSX.Fuselage_flags.STL_fus;
rotate_thetax = 180; % rotates the x axis
rotate_thetay = 270; % rotates the y axis
% [original, movelist, z_slices_tmp] = STL_slicer(STL_fus,slice_height,rotate_thetax,rotate_thetay,OUTPUT_read_XLSX);
[original, movelist, z_slices_tmp] = STL_slicer_v3(STL_fus,slice_height,rotate_thetax,rotate_thetay,OUTPUT_read_XLSX);
% plot_slices(movelist,z_slices_tmp, 0)
% pause

model = stl2matlab(OUTPUT_read_XLSX.Fuselage_flags.STL_fus);
ScaleFactor = OUTPUT_read_XLSX.AC_Data_flags.SF; % CAD designs are in mm, and changing to m
SF_CAD = OUTPUT_read_XLSX.Fuselage_flags.SF_CAD;

plotData{1} = (model{1}*ScaleFactor).*SF_CAD;
plotData{2} = (model{2}*ScaleFactor).*SF_CAD;
plotData{3} = (model{3}*ScaleFactor).*SF_CAD;

% plots to check if fuselages has been ortated accordingly
PLOT_Slices_fuselage_STL = OUTPUT_read_XLSX.Fuselage_flags.PLOT_Slices_fuselage_STL;

% %% Asigns the values obtained from STL file
% % Eliminates empty Cells
% check_empy = (~cellfun('isempty',movelist));
% j=1;
% for i = 1:length(movelist)
%     if check_empy(i) == 1;
%         % reconstruct the scanned data removing empty cells
%         movelist_mod{j} = movelist{i};
%         z_slices(j) = z_slices_tmp(i);
%     j = j+1;
%     end
% end


%% Asigns the values obtained from STL file
% Eliminates empty Cells
check_empy = (~cellfun('isempty',movelist));
j=1;
for i = 1:length(movelist)
    if check_empy(i) == 0
        if i==1
            movelist_mod{i} = movelist{2}*0.5;
        else
            check_low = movelist{i-1}(1,1);
            check_high = movelist{i+i}(1,1);
            factor_b = check_low/check_high;
            movelist_mod{i} = movelist{i+1}*factor_b;
            movelist{1} = movelist_mod{i};
            % z_slices(i) = z_slices_tmp(i);
        end              
    else
        movelist_mod{i} = movelist{i};
    end
end

% includes the origin to the z_slices vector
z_slices = [0 z_slices_tmp];

% movelist = (~cellfun('isempty',movelist));
% check_movelist_mod = (~cellfun('isempty',movelist_mod));
movelist = movelist_mod;

for i=1:length(movelist_mod)
    % construct the first point using the first value of the slices and    % reduce it
    % Eliminates NaN
    Y = movelist_mod{i}(:,2)*SF_CAD;
    Y = Y';
    Y = Y(~isnan(Y))';
    Y_stl{i+1}(1,:) = Y;
    % Eliminates NaN
    Z = movelist_mod{i}(:,1)*SF_CAD;
    Z = Z';
    Z = Z(~isnan(Z))';
    Z_stl{i+1}(1,:) = Z;
    J = length(movelist_mod{i});
    for j=1:J
        X_stl{i+1}(1,j) = z_slices(i+1)*SF_CAD;
    end
end

% Stores the origonal not shifted
X_stl_original = X_stl;
Y_stl_original = Y_stl;
Z_stl_original = Z_stl;

% Calculates offset value
minY = min(Y_stl{3});
maxY = max(Y_stl{3});
minZ = min(Z_stl{3});
maxZ = max(Z_stl{3});
Delta_Y_stl = maxY - (maxY - minY)/2;
Delta_Z_stl = maxZ - (maxZ - minZ)/2;

% Shifts Y and Z data with the offset valu such thta the nose is the 0,0,0
for i=1:length(X_stl)
    Y_stl{i} = Y_stl{i} - Delta_Y_stl;
    Z_stl{i} = (Z_stl{i} - Delta_Z_stl);
end

Body_Geo.Delta_Z_stl = Delta_Y_stl;
Body_Geo.Delta_Y_stl = Delta_Z_stl;

% Augments the vectors by Assigning values for the first element using the values of the second element
red = 0.5;
Y_stl{2}(1,:) = Y_stl{3}(1,:)*red;
Z_stl{2}(1,:) = Z_stl{3}(1,:)*red;
J = length(Z_stl{3}(1,:));
for j=1:J
    X_stl{2}(1,j) = z_slices(2)*SF_CAD;
end

red = 0.5;
Y_stl{1}(1,:) = Y_stl{2}(1,:)*red;
Z_stl{1}(1,:) = Z_stl{2}(1,:)*red;
J = length(Z_stl{2}(1,:));
for j=1:J
    X_stl{1}(1,j) = z_slices(1)*SF_CAD;
end

% Augments the vectors by Assigning values for the last element using the values of the last  element
red = 1;
END_value = length(z_slices);
Y_stl{END_value}(1,:) = Y_stl{end}(1,:)*red;
Z_stl{END_value}(1,:) = Z_stl{end}(1,:)*red;
J = length(Z_stl{END_value}(1,:));
for j=1:J
    X_stl{END_value}(1,j) = z_slices(end)*SF_CAD;
end

% % Augments the vectors by Assigning values for the first element using the values of the second element
% red = 0.5;
% Y_stl{1}(1,:) = Y_stl{2}(1,:)*red;
% Z_stl{1}(1,:) = Z_stl{2}(1,:)*red;
% J = length(Z_stl{1}(1,:));
% for j=1:J
%     X_stl{1}(1,j) = z_slices(1);
% end

% % Calculates offset value
% minY = min(Y_stl{3});
% maxY = max(Y_stl{3});
% minZ = min(Z_stl{3});
% maxZ = max(Z_stl{3});
% Delta_Y_stl = maxY - (maxY - minY)/2;
% Delta_Z_stl = maxZ - (maxZ - minZ)/2;
% % Shifts Y and Z data with the offset valu such thta the nose is the 0,0,0
% for i=1:length(X_stl)
%     Y_stl{i} = Y_stl{i} - Delta_Y_stl;
%     Z_stl{i} = (Z_stl{i} - Delta_Z_stl);
% end

% % Plots only if want to check geometry
% plots = OUTPUT_read_XLSX.PLOT_flags.plot;
% if plots(5)==1
%     if PLOT_Slices_fuselage_STL == 1
%         Fig = Fig + 1;
%         figure(Fig)
%         subplot(1,3,1)
%         plot_slices(movelist,z_slices,0);
%         grid on
%         title('Fuselage rotated around the y-axis to be sliced')
%         xlabel('Fuselage z-axis')
%         ylabel('Fuselage y-axis')
%         zlabel('Fuselage x-axis')
% 
%         subplot(1,3,2)
%         for i=3:length(z_slices)-1
%             plot(Y_stl_original{i},Z_stl_original{i})
%             hold on
%         end
%         title('Fuselage y-axis vs z-axis - Not Shifted (original data)')
%         xlabel('Fuselage y-axis')
%         ylabel('Fuselage z-axis')
%         grid on
%         hold off
% 
%         subplot(1,3,3)
%         for i=1:length(z_slices)
%             plot(Y_stl{i},Z_stl{i})
%             hold on
%         end
%         title('Fuselage y-axis vs z-axis - Shifted to origin (@ the nose)')
%         xlabel('Fuselage y-axis')
%         ylabel('Fuselage z-axis')
%         grid on
%         hold off
%         %         prompt = sprintf('CHECK IF THE FUSELAGE Y AND Z AXIS ARE DRAWN CORRECTY, IF SO RUNTHE SIMULATION AGAIN SETTING "PLOT_Slices_fuselage_STL=0" IN SHEET 6 OF EXCEL');
%     end
% end

% % Calculo de los perímetros de cada sección
% for i=1:(n_sections_body)
%     Per(i).y = [Data_P(i).y; -Data_P(i).y;Data_P(i).y(1)];
%     Per(i).z = [Data_P(i).z; fliplr(Data_P(i).z')';Data_P(i).z(1)];
% end

% Calculo de los perímetros de cada sección
for i=1:length(X_stl)
    Per(i).y = Y_stl{i};
    Per(i).z = Z_stl{i};
    x_max(i) = max(X_stl{i});
    y_max(i) = max(Y_stl{i});
    z_max(i) = max(Z_stl{i});
end

% plotData{1}  = [x(:,end:(-1):1),x];
% plotData{2}  = [-y(:,end:(-1):1), y];
% plotData{3} = [z(:,end:(-1):1), z];

% figure(3)
% mesh(plotData{1},plotData{2},plotData{3});
% pause

% for i=1:(n_sections_body-2)
for i=1:length(z_slices)
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
alpha               = -atan(Delta_h./Delta_Zcent);
Body_Geo.alpha_x    = [alpha,alpha(end)];
% Calculo del volumen (aproximación)
l=1;
% for i=3:n_sections_body-1
for i=2:length(z_slices)
    z_slices_surface(i-1) = z_slices(i);
    Delta_h = x_max(i) - x_max(i-1);
    Surface(l) = Delta_h*(Perimeter(i) + Perimeter(i-1))/2;
    l = l+1;
end

% Surface total geometría
Surf_TOT = sum(Surface);
Body_Geo.Surf_TOT = Surf_TOT;

% Calculo de las áreas de cada sección
% for i=1:n_sections_body
for i=1:length(z_slices)
    %     Area_body(i) = polyarea(Data_P(i).y,Data_P(i).z);
    Area_body(i) = polyarea(Y_stl{i},Z_stl{i});
    %     dSdX(i) = diff(Area_body(i))./diff(x_max(i));
end

% Body_Geo.Area_body = Area_body;
Body_Geo.Area_body = Area_body;

% X-positions of fuselage
Body_Geo.x_Area_body = x_max;

% Calculo de los perímetros de cada sección
% for i=1:(n_sections_body)
for i=1:length(z_slices)
    Height(i) = max(Per(i).z) - min(Per(i).z);
    Bisectriz_z(i) = Height(i)/2;
    Centroid_z(i) =  max(Per(i).z) - Bisectriz_z(i);
end

% if plots(5)==1
%     if PLOT_Slices_fuselage_STL == 1
%         Fig = Fig + 1;
%         figure(Fig)
%         subplot(2,2,1)
%         plot(z_slices,Area)
%         grid on
%         title('Fuselage Area persections')
%         xlabel('Fuselage x-axis sections')
%         ylabel('Fuselage Area per section')
% 
%         subplot(2,2,2)
%         plot(z_slices,Area_body)
%         grid on
%         title('Fuselage Area persections (method 2)')
%         xlabel('Fuselage x-axis sections')
%         ylabel('Fuselage Area per section')
% 
%         subplot(2,2,3)
%         plot(z_slices,Perimeter)
%         grid on
%         title('Fuselage Area persections (method 2)')
%         xlabel('Fuselage x-axis sections')
%         ylabel('Fuselage Area per section')
% 
%         subplot(2,2,4)
%         plot(z_slices_surface,Surface)
%         grid on
%         title('Surface Area per sections (pairs of sections)')
%         xlabel('Fuselage x-axis sections')
%         ylabel('Fuselage Surface per sections')
% 
%         Fig = Fig + 1;
%         figure(Fig)
%         subplot(2,3,1)
%         plot(z_slices,X_cen)
%         grid on
%         title('Fuselage Centroid of each Area persections')
%         xlabel('Fuselage x-axis sections')
%         ylabel('Fuselage X-Centroid ')
% 
%         subplot(2,3,2)
%         plot(z_slices,Y_cen)
%         grid on
%         title('Fuselage Centroid of each Area persections')
%         xlabel('Fuselage x-axis sections')
%         ylabel('Fuselage Y-Centroid')
% 
%         subplot(2,3,3)
%         plot(z_slices,Z_cen)
%         grid on
%         title('Fuselage Centroid of each Area persections')
%         xlabel('Fuselage x-axis sections')
%         ylabel('Fuselage Z-Centroid')
% 
%         subplot(2,3,4)
%         plot(z_slices,Height)
%         grid on
%         title('Fuselage Heigth of each Area per sections')
%         xlabel('Fuselage x-axis sections')
%         ylabel('Fuselage Height ')
% 
%         subplot(2,3,5)
%         plot(z_slices,Bisectriz_z)
%         grid on
%         title('Fuselage Bisectric of each Area persections')
%         xlabel('Fuselage x-axis sections')
%         ylabel('Fuselage Bisectriz')
% 
%         subplot(2,3,6)
%         plot(z_slices,Centroid_z)
%         grid on
%         title('Fuselage Centroid of each Area persections')
%         xlabel('Fuselage x-axis sections')
%         ylabel('Fuselage Centroid')
%     end
% end

Angle_fus = gradient(Centroid_z,x_max)*R2D;
x_interp = linspace(x_max(1),x_max(end),200);
Centroid_z_interp = interp1(x_max,Centroid_z,x_interp,'spline');
Angle_fus_interp = gradient(Centroid_z_interp,x_interp)*R2D;

Body_Geo.Angle_fus = Angle_fus;
Body_Geo.Angle_fus_interp = Angle_fus_interp;
Body_Geo.x_interp = x_interp;

% Calculo del ángulo de cada sección
% for i=2:n_sections_body
for i=2:length(z_slices)
    Delta_h(i-1) = x_max(i) - x_max(i-1);
    %     Angle_fus(i-1) = atan(Centroid_z(i)/Delta_h(i-1));
end

Body_Geo.Height = Height;
Body_Geo.Bisectriz_z = Bisectriz_z;
Body_Geo.Centroid_z = Centroid_z;
Body_Geo.Angle_fus = Angle_fus;

% Calculo del volumen (aproximación)
l=1;
% for i=2:n_sections_body
for i=2:length(z_slices)
    z_slices_surface(i-1) = z_slices(i);
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
Area_top  = polyarea(x_max,y_max);
% Calculo del Area proyectada en plano x-z
Area_side = polyarea(x_max,z_max);

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

i_SfrontMax = 2; %  does not include the firs elements
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
Body_Geo.H_x = 2*z_max;
Body_Geo.h_fus   = max(Body_Geo.H_x);
% Diameter
Body_Geo.D_x = (Body_Geo.H_x + Body_Geo.W_x)/2;
Body_Geo.d_fus   = max(Body_Geo.D_x);

Body_Geo.x_max   = x_max;
Body_Geo.y_max   = y_max;
Body_Geo.z_max   = z_max;

%% PROPIEDADES AERODINÁMICAS
% Cálculo de la pendiente de sustentación del fuselaje
Fineness_Ratio = Body_Geo.l_fus/Body_Geo.w_Area_b_max;
%digitaliazacion figura PAMADI CAP3 Fig 3.6
x_f_k2_k1 = [4.,5.,6.,8.,10.,12.,14.,16.,18.,20.];
y_f_k2_k1 = [.77,.825,.865,.91,.94,.955,.965,.97,.973,.975];
f_k2_k1  = interp1(x_f_k2_k1,y_f_k2_k1,Fineness_Ratio,'spline');
% f_k2_k1         = k2_k1_calc(Body_Geo.l_fus, Geo_tier.d_fus);      % Fuselage apparent mass coefficient. Pamadi, Figure 3.6
Body_Geo.CLa_fus    = 2*f_k2_k1*(Body_Geo.S_front/Sref);         % Pendiente de sustentación del morro aislado
test1=0;
% if plots(5)==1
%     if PLOT_Slices_fuselage_STL == 1
%         Fig = Fig + 1;
%         figure(Fig)
%         subplot(1,3,1)
%         plot(z_slices,Angle_fus)
%         hold on
%         grid on
%         plot(x_interp,Angle_fus_interp)
%         legend('angle','angle interpolated')
%         title('Fuselage Angle persections')
%         xlabel('Fuselage x-axis sections')
%         ylabel('Fuselage Angle per section (degs)')
% 
%         subplot(1,3,2)
%         plot(z_slices_surface,Vol)
%         grid on
%         title('Fuselage Volume per sections')
%         xlabel('Fuselage x-axis sections')
%         ylabel('Fuselage Volume per section (m^3)')
% 
%         subplot(1,3,3)
%         plot(z_slices(2:end),dSdX(2:end))
%         grid on
%         title('Fuselage dS/dx per sections')
%         xlabel('Fuselage x-axis sections')
%         ylabel('Fuselage dS/dx per section')
%     end
% end
% Fig = Fig + 1;
end