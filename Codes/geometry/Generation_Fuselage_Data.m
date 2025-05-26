function [Body_Geo,meshData] = Generation_Fuselage_Data(Geo_tier,OUTPUT_read_XLSX,XFLR5_DATA,XFLR5_file,STL_file,conversion_UNITS)  % Defines Fuselage DATA

R2D = conversion_UNITS.R2D;

% Retieving data
scale = OUTPUT_read_XLSX.AC_Data_flags.SF;
CASE_fuse = OUTPUT_read_XLSX.AC_Data_flags.CASE_fuse;
ESCALADO = OUTPUT_read_XLSX.AC_Data_flags.ESCALADO;
% XFLR5_DATA = OUTPUT_read_XLSX.Fuselage_flags.XFLR5_DATA;
% XFLR5_file = OUTPUT_read_XLSX.Fuselage_flags.XFLR5_file;
% STL_file = OUTPUT_read_XLSX.Fuselage_flags.STL_file;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates Geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% XFLR5_file = strcat('fuse_pepi_desplazado.txt');
% XFLR5_file = strcat('Fuselage_DATA_scaled.txt');
% STL_file = strcat('Fus_ProVant_v4.0.stl');

%--------------------------- CALCULATES GEOMETRY DATA ---------------------------------
% Data to read XFLR5 files
nSections = XFLR5_DATA.nSections;
nPoints = XFLR5_DATA.nPoints;
lecture_Geo = XFLR5_DATA.lecture_Geo;
STL_PLOT = XFLR5_DATA.STL_PLOT;

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

switch CASE_fuse
    case 1 % STL of Fuselaje
        Ratio_SCALE_l_fus = Ratio_SCALE_tmp.Ratio_SCALE_l_fus;
        Ratio_SCALE_w_fus = Ratio_SCALE_tmp.Ratio_SCALE_w_fus;
        Ratio_SCALE_h_fus = Ratio_SCALE_tmp.Ratio_SCALE_h_fus;
        
        file_stl = 'ProVant1.stl';
        triangles = read_binary_stl_file(file_stl);
        original = triangles;
        triangles = rotate_stl(triangles,'y',90);
        triangles = rotate_stl(triangles,'z',180);
        slice_height = 1;
        
        % Scaling Factor
        ScalingFactor = 1/1000;
        [movelist, yz_slices] = slice_stl_create_path(triangles, slice_height);
        
        % Scale from mm to meters
        yz_slices = ScalingFactor*yz_slices;
        
        size_slices = size(yz_slices);
        n_SLICES = size_slices(2);
        for i = 1:(n_SLICES)
            x(i) = (yz_slices(i) - yz_slices(1));
        end
        
        size_movelist = size(movelist);
        n_SETS = size_movelist(2);
        
        % Estimates the offset value using the second slice of the set
        % Remove NaN
        Y_off = ScalingFactor*movelist{2}(:,2);
        Y_off= Y_off(~isnan(Y_off));
        % Remove NaN
        Z_off = ScalingFactor*movelist{2}(:,1);
        Z_off = Z_off(~isnan(Z_off));
        [geom_offset,iner_offset,cpmo_offset] = polygeom(Y_off,Z_off);
        %   GEOM = [ area   X_cen  Y_cen  perimeter ]
        %   INER = [ Ixx    Iyy    Ixy    Iuu    Ivv    Iuv ]
        %     u,v are centroidal axes parallel to x,y axes.
        %   CPMO = [ I1     ang1   I2     ang2   J ]
        %     I1,I2 are centroidal principal moments about axes
        %         at angles ang1,ang2.
        %     ang1 and ang2 are in radians.
        %     J is centroidal polar moment.  J = I1 + I2 = Iuu + Ivv
        % shift the values to the origin
        Y_offset = geom_offset(2);
        Z_offset = geom_offset(3);
        
        movelist{1} = [0 0];
        % Scaling from mm to meters
        Y{1} = ScalingFactor*movelist{1}(:,2);
        Z{1} = ScalingFactor*movelist{1}(:,1);
        X{1} = x(1);
        Max_size = 0;
        for i = 2:(n_SETS)
            % Remove NaN & offset
            Y{i} = ScalingFactor*movelist{i}(:,2);
            Y{i}= Y{i}(~isnan(Y{i}))- Y_offset;
            size_Y = size(Y{i});
            % Remove NaN
            Z{i} = ScalingFactor*movelist{i}(:,1);
            Z{i}= Z{i}(~isnan(Z{i}))- Z_offset;
            
            X{i} = x(i)*ones(size_Y);
            [geom{i},iner{i},cpmo{i}] = polygeom(Y{i},Z{i});
            %   GEOM = [ area   X_cen  Y_cen  perimeter ]
            %   INER = [ Ixx    Iyy    Ixy    Iuu    Ivv    Iuv ]
            %     u,v are centroidal axes parallel to x,y axes.
            %   CPMO = [ I1     ang1   I2     ang2   J ]
            %     I1,I2 are centroidal principal moments about axes
            %         at angles ang1,ang2.
            %     ang1 and ang2 are in radians.
            %     J is centroidal polar moment.  J = I1 + I2 = Iuu + Ivv
            Area(i) = geom{i}(1);
            X_cen(i) = x(i);
            Y_cen(i) = geom{i}(2);
            Z_cen(i) = geom{i}(3);
            Perimeter(i) = geom{i}(4);
            
            x_max(i) = max(X{i})*Ratio_SCALE_l_fus;
            y_max(i) = max(Y{i})*Ratio_SCALE_w_fus;
            z_max(i) = max(Z{i})*Ratio_SCALE_h_fus;

        end
        
        Body_Geo.S_x    = Area;
        Body_Geo.x      = x_max;
        
        Delta_h             = Body_Geo.x(2:end) - Body_Geo.x(1:(end-1));

        Delta_Zcent         = Z_cen(2:end)-Z_cen(1:(end-1));
        alpha               = -atan(Delta_h./Delta_Zcent);
        Body_Geo.alpha_x    = [alpha,alpha(end)];

        l=1;
        for i=3:n_SETS-1
            Delta_h = x_max(i) - x_max(i-1);
            Surface(l) = Delta_h*(Perimeter(i-1) + Perimeter(i-2))/2;
            l = l+1;
        end
        
        % Surface total geometr?a
        Surf_TOT = sum(Surface);
        Body_Geo.Surf_TOT = Surf_TOT;
        
        % Calculo de las áreas de cada sección
        Area_body = Area;
        Body_Geo.Area_body = Area_body;
        
        % X-positions of fuselage
        Body_Geo.x_Area_body = x_max;
        
        % Calculo de los perímetros de cada sección
        for i=1:(n_SETS)
            Height(i) = max(Z{i}) - min(Z{i});
            Bisectriz_z(i) = Height(i)/2;
            Centroid_z(i) =  max(Z{i}) - Bisectriz_z(i);
        end
        
        % Calculo de los perímetros de cada sección
        Angle_fus = gradient(Z_cen,x_max)*R2D;
        
        Body_Geo.Angle_fus = Angle_fus;
        Body_Geo.Angle_fus_interp = Angle_fus;
        Body_Geo.x_interp = x_max;
        
        % Calculo del ángulo de cada sección
        for i=2:n_SETS
            Delta_h(i-1) = x_max(i) - x_max(i-1);
            %     Angle_fus(i-1) = atan(Centroid_z(i)/Delta_h(i-1));
        end
        
        Body_Geo.Height = Height;
        Body_Geo.Bisectriz_z = Bisectriz_z;
        Body_Geo.Centroid_z = Z_cen;
        
        % Calculo del volumen (aproximación)
        l=1;
        for i=2:n_SETS
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
        Body_Geo.H_x = 2*z_max;
        Body_Geo.h_fus   = max(Body_Geo.H_x);
        % Diameter
        Body_Geo.D_x = (Body_Geo.H_x + Body_Geo.W_x)/2;
        Body_Geo.d_fus   = max(Body_Geo.D_x);
        
        Body_Geo.x_max   = x_max;
        Body_Geo.y_max   = y_max;
        Body_Geo.z_max   = z_max;
        
        % PROPIEDADES AERODINÁMICAS
        % Cálculo de la pendiente de sustentación del fuselaje
%         f_k2_k1         = k2_k1_calc(Body_Geo.l_fus, Geo_tier.d_fus);      % Fuselage apparent mass coefficient. Pamadi, Figure 3.6
%         Body_Geo.CLa_fus    = 2*f_k2_k1*(Body_Geo.S_front/Sref);         % Pendiente de sustentación del morro aislado

    case 2
        % sTL of complete aircraft
        [tri, fileform, A, S] = stlread('ProVant.stl')
        
    case 3 % Read File from XFLR5
        
        Ratio_SCALE_l_fus = Ratio_SCALE_tmp.Ratio_SCALE_l_fus;
        Ratio_SCALE_w_fus = Ratio_SCALE_tmp.Ratio_SCALE_w_fus;
        Ratio_SCALE_h_fus = Ratio_SCALE_tmp.Ratio_SCALE_h_fus;

        n_ele_geo = nPoints; % elementos que se han de leer
        n_sections_body = nSections; % Elemtos del fuselaje

        for i=1:n_sections_body
            lecture_Geo = lecture_Geo;
            % ExtracciÃ³n de datos VLM
            file_id=fopen(XFLR5_file);
            aux = textscan(file_id,'%n%n%n','Headerlines',lecture_Geo);
            %      size(aux{i})
            %      pause
            Data_P(i).raw_data= aux;
            Data_P(i).label = XFLR5_file;
            Data_P(i).x = aux{1}*Ratio_SCALE_l_fus;
            Data_P(i).y = aux{2}*Ratio_SCALE_w_fus;
            Data_P(i).z = aux{3}*Ratio_SCALE_h_fus;
            x(:,i)      = aux{1}*Ratio_SCALE_l_fus;
            y(:,i)      = aux{2}*Ratio_SCALE_w_fus;
            z(:,i)      = aux{3}*Ratio_SCALE_h_fus;
            x_max(i) = max(aux{1}*Ratio_SCALE_l_fus);
            y_max(i) = max(aux{2}*Ratio_SCALE_w_fus);
            z_max(i) = max(aux{3}*Ratio_SCALE_h_fus);
            fclose(file_id);
            lecture_Geo = lecture_Geo + n_ele_geo  + 2;
        end 
        
        % Calculo de los perï¿½metros de cada secciï¿½n
        for i=1:(n_sections_body)
            Per(i).y = [Data_P(i).y; -Data_P(i).y;Data_P(i).y(1)];
            Per(i).z = [Data_P(i).z; fliplr(Data_P(i).z')';Data_P(i).z(1)];
        end
        
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
        alpha               = -atan(Delta_h./Delta_Zcent);
        Body_Geo.alpha_x    = [alpha,alpha(end)];
        
        % Calculo del volumen (aproximaciÃ³n)
        l=1;
        for i=3:n_sections_body-1
            Delta_h = x_max(i) - x_max(i-1);
            Surface(l) = Delta_h*(Perimeter(i-1) + Perimeter(i-2))/2;
            l = l+1;
        end
        
        % Surface total geometrï¿½a
        Surf_TOT = sum(Surface);
        Body_Geo.Surf_TOT = Surf_TOT;
        
        % Calculo de las ï¿½reas de cada secciï¿½n
        for i=1:n_sections_body
            Area_body(i) = 2*polyarea(Data_P(i).y,Data_P(i).z);
            %     dSdX(i) = diff(Area_body(i))/diff(x_max(i));
        end
        
        % Body_Geo.Area_body = Area_body;
        Body_Geo.Area_body = Area_body;
        
        % X-positions of fuselage
        Body_Geo.x_Area_body = x_max;
        
        % Calculo de los perï¿½metros de cada secciï¿½n
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
        f_k2_k1         = k2_k1_calc(Body_Geo.l_fus, Geo_tier.d_fus);      % Fuselage apparent mass coefficient. Pamadi, Figure 3.6
        Body_Geo.CLa_fus    = 2*f_k2_k1*(Body_Geo.S_front/Geo_tier.S_ref);         % Pendiente de sustentación del morro aislado
        
%         Body_Geo
%         pause
        plotData{1}  = [x(:,end:(-1):1),x];
        plotData{2}  = [-y(:,end:(-1):1), y];
        plotData{3} = [z(:,end:(-1):1), z];

        DATA_FUS{1} = x;
        DATA_FUS{2} = y;
        DATA_FUS{3} = z;
     
end


% Variable that define the scaling factor
Ratio_SCALE.Ratio_SCALE_l_fus = Ratio_SCALE_l_fus;
Ratio_SCALE.Ratio_SCALE_w_fus = Ratio_SCALE_w_fus;
Ratio_SCALE.Ratio_SCALE_h_fus = Ratio_SCALE_h_fus;

% [Body_Geo_XFLR5,meshData_XFLR5] = Calc_Body_Geometry_Feb2021(XFLR5_file,STL_file,Geo_tier, nPoints, nSections,conversion_UNITS,lecture_Geo,Ratio_SCALE_tmp,scale,plotData,DATA_FUS);
[Body_Geo,meshData] = Calc_Body_Geometry_Dec2018(XFLR5_file,STL_file,Geo_tier, nPoints, nSections,conversion_UNITS,lecture_Geo,Ratio_SCALE,scale);

% Scale from fuselage in XFLR5 to CAD
l_fus1 = Body_Geo_XFLR5.l_fus; % Length of original XFLR5
w_fus1 = Body_Geo_XFLR5.w_fus; % width
h_fus1 = Body_Geo_XFLR5.h_fus; % height

% Desired Lenth from input data
l_fus2 = Geo_tier.l_fus; % Desired Length
w_fus2 = Geo_tier.w_fus; % width
h_fus2 = Geo_tier.h_fus; % height
% Determine the type of scaling
% ESCALADO = 0; % no scaling Ratio =1
% ESCALADO = 1; % scaling acording to desired for each 3 axix
% ESCALADO = 2; % Scaling for emergentia, maintaing the same ratio all 3 axis
if ESCALADO == 0
    % Maintaining the ratio
    Ratio_SCALE_l_fus = 1;
    Ratio_SCALE_w_fus = 1;
    Ratio_SCALE_h_fus = 1;
elseif ESCALADO == 1
    % Applies scaling factor according to desired geometry
    Ratio_SCALE_l_fus = l_fus2/l_fus1;
    Ratio_SCALE_w_fus = w_fus2/w_fus1;
    Ratio_SCALE_h_fus = h_fus2/h_fus1;
elseif ESCALADO == 2
    % Maintaining the ratio
    Ratio_SCALE_l_fus = l_fus2/l_fus1;
    Ratio_SCALE_w_fus = l_fus2/l_fus1;
    Ratio_SCALE_h_fus = l_fus2/l_fus1;
end

% Variable that define the scaling factor
Ratio_SCALE.Ratio_SCALE_l_fus = Ratio_SCALE_l_fus;
Ratio_SCALE.Ratio_SCALE_w_fus = Ratio_SCALE_w_fus;
Ratio_SCALE.Ratio_SCALE_h_fus = Ratio_SCALE_h_fus;

% Being used only if STL model used to determine geometry of fuselage
% Nth = size(plotData{1});
% STLnSections = Nth(1);% number of sectins (x coordinate) (eliminated 1 for convergence)
% STLnPoints = Nth(2); % number of elements per section (y and z coordinates)
% Determines fuselage information
[Body_Geo,meshData] = Calc_Body_Geometry_Dec2018(XFLR5_file,STL_file,Geo_tier, nPoints, nSections,conversion_UNITS,lecture_Geo,Ratio_SCALE,scale);
% Identifies the fuselage mode, between data from XFLR5 and STL files (CAD)
if STL_PLOT == 1
    switch CASE_fuse
        case 1 % STL of Fuselaje
            model = stl2matlab('ProVant.stl');
            ScaleFactor = scale*1/1000; % CAD designs are in mm, and changing to m
            plotData{1} = model{1}*ScaleFactor + Geo_tier.x_offset_CAD;
            plotData{2} = model{2}*ScaleFactor;
            plotData{3} = model{3}*ScaleFactor + Geo_tier.z_offset_CAD;
            meshData = plotData;
        case 2
            % sTL of complete aircraft
            model = stl2matlab('VTOL_Single.stl');
            model = stl2matlab('VTP_WIG.stl');
            model = stl2matlab('VTP_WIG.stl');            
            ScaleFactor = scale; % CAD designs are in mm, and changing to m
            plotData{1} = model{1}*ScaleFactor + Geo_tier.x_offset_CAD;
            plotData{2} = model{2}*ScaleFactor;
            plotData{3} = model{3}*ScaleFactor + Geo_tier.z_offset_CAD;
            meshData = plotData;

            model = stl2matlab('HTP_WIG.stl');            
            ScaleFactor = scale; % CAD designs are in mm, and changing to m
            plotData{1} = model{1}*ScaleFactor + Geo_tier.x_offset_CAD;
            plotData{2} = model{2}*ScaleFactor;
            plotData{3} = model{3}*ScaleFactor + Geo_tier.z_offset_CAD;
            meshData = plotData;
        case 3
            meshData = meshData;    
    end    
end