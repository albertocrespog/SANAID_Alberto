function Fig = plot_GEOMETRY_2022(Geo_tier,Plot_Options,Body_Geo,meshData,Prop_data,Fig,AC_CONFIGURATION,COLOR_scheme,case_AC,Performance,OUTPUT_read_XLSX,filenameS)

% Aircraft type
% AC_type = 1 - flying wing
% AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
% AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
% AC_type = 4 - 2 surface: wing + V-tail
% AC_type = 5 - 3 surface: cannard + wing + V-tail
% Engine location
% Engine_loc = 1 - under wings
% Engine_loc = 2 - fuselage
% Engine_loc = 3 - wingtips
% Engine Configuration
% Engine_conf = 1 - pusher prop
% Engine_conf = 2 - puller prop
% Engine_conf = 3 - turbine

% Defines the Target for the plots
fname = filenameS.filename_Plots;

% Stores Aircraft Configuration
% AC_CONFIGURATION.Control_surface = Control_surface;
AC_type = AC_CONFIGURATION.AC_type;
Engine_loc = AC_CONFIGURATION.Engine_loc;
Engine_conf = AC_CONFIGURATION.Engine_conf;

W1 = AC_CONFIGURATION.W1;
HTP = AC_CONFIGURATION.HTP;
VTP = AC_CONFIGURATION.VTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
twin_VTP = AC_CONFIGURATION.twin_VTP;

% Plotting options
FS = Plot_Options.FS;
Video_3D = Plot_Options.Video_3D;
SAVE_FIGS = Plot_Options.SAVE_FIGS;
L_fus = Geo_tier.l_fus;
n_eng = OUTPUT_read_XLSX.Propulsive_flags.propul(2);

filename_Plots = filenameS.filename_Plots;
ESCALADO = OUTPUT_read_XLSX.AC_Data_flags.ESCALADO;
[XFLR5_DATA,XFLR5_file,STL_files] = Extraction_fuselage_data(OUTPUT_read_XLSX);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates Geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% XFLR5_file = strcat('fuse_pepi_desplazado.txt');
% XFLR5_file = strcat('Fuselage_DATA_scaled.txt');
% STL_file = strcat('Fus_ProVant_v4.0.stl');
SF = OUTPUT_read_XLSX.AC_Data_flags.SF;
CASE_fuse = OUTPUT_read_XLSX.AC_Data_flags.CASE_fuse;

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
% Reads to extract length of fuselag e 
% assumes Ratio_SCALE = 1
Ratio_SCALE_l_fus_tmp = 1;
Ratio_SCALE_w_fus_tmp = 1;
Ratio_SCALE_h_fus_tmp = 1;
Ratio_SCALE_tmp.Ratio_SCALE_l_fus = Ratio_SCALE_l_fus_tmp; 
Ratio_SCALE_tmp.Ratio_SCALE_w_fus = Ratio_SCALE_w_fus_tmp;
Ratio_SCALE_tmp.Ratio_SCALE_h_fus = Ratio_SCALE_h_fus_tmp;
% [Body_Geo_XFLR5,meshData_XFLR5] = Calc_Body_Geometry_Dec2018(XFLR5_file,STL_files,Geo_tier, nPoints, nSections,conversion_UNITS,lecture_Geo,Ratio_SCALE_tmp,SF,OUTPUT_read_XLSX);

STL_fus = OUTPUT_read_XLSX.Fuselage_flags.STL_fus;
SF_CAD = OUTPUT_read_XLSX.Fuselage_flags.SF_CAD;
% Reads data for fuselage
[VertexData,FVCD,isBinary] = stlread(STL_fus);
% Scale from fuselage read from STL
l_fus1 = (max(VertexData.Points(:,1)) - min(VertexData.Points(:,1)))*SF_CAD;
w_fus1 = (max(VertexData.Points(:,2)) - min(VertexData.Points(:,2)))*SF_CAD;
h_fus1 = (max(VertexData.Points(:,3)) - min(VertexData.Points(:,3)))*SF_CAD;

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

[Fig,Body_Geo,plotData] = Generate_Plots_Body_Geometry(XFLR5_file,STL_files,Geo_tier,nPoints, nSections,conversion_UNITS,lecture_Geo,Ratio_SCALE,SF,OUTPUT_read_XLSX,Fig);

% Fig = Fig +1;
% figure(Fig)
PLOTTING_UAV = plot_Models_ITER(Geo_tier,Body_Geo,meshData,Prop_data,AC_CONFIGURATION,OUTPUT_read_XLSX);

%% plotting fuselage
meshData = PLOTTING_UAV.meshData; % FUSELAJE

%% Defines colors RGB
color_fus = COLOR_scheme.color_fus;
color_w1 = COLOR_scheme.color_w1;
color_HTP = COLOR_scheme.color_HTP;
color_vee = COLOR_scheme.color_vee;
color_vee2 = COLOR_scheme.color_vee2;
color_can = COLOR_scheme.color_can;
color_prop = COLOR_scheme.color_prop;
color_vtp = COLOR_scheme.color_vtp;
color_eng = COLOR_scheme.color_eng;
color_nac = COLOR_scheme.color_nac;
color_ac = COLOR_scheme.color_ac;

Delta_Z_stl = Body_Geo.Delta_Z_stl;
Delta_Y_stl = Body_Geo.Delta_Y_stl;

% % Calculates the Aircraft Geometry
% model = stl2matlab(OUTPUT_read_XLSX.Fuselage_flags.STL_ac);
% ScaleFactor = OUTPUT_read_XLSX.AC_Data_flags.SF; % CAD designs are in mm, and changing to m
% SF_CAD = OUTPUT_read_XLSX.Fuselage_flags.SF_CAD;
% [points_ac,triangles,tri norms] = import_stl_fast(OUTPUT_read_XLSX.Fuselage_flags.STL_ac,1);
% 
% max_x_ac = max(points_ac(:,1));
% min_x_ac = min(points_ac(:,1));
% max_y_ac = max(points_ac(:,2));
% min_y_ac = min(points_ac(:,2));
% max_z_ac = max(points_ac(:,3));
% min_z_ac = min(points_ac(:,3));
% 
% Delta_X_stl = points_ac(1,1);
% Delta_Y_stl = points_ac(1,2);
% Delta_Z_stl = points_ac(1,3);
% 
% % Offsets
% % plotData{1} = (model{1}*ScaleFactor + Geo_tier.x_offset_CAD - min_x_ac).*SF_CAD;
% % plotData{2} = (model{2}*ScaleFactor - (max_y_ac+min_y_ac)/2).*SF_CAD;
% % plotData{3} = (model{3}*ScaleFactor + Geo_tier.z_offset_CAD - Delta_Z_stl).*SF_CAD;
% % meshData_ac = plotData;
%
% plotData{1} = (model{1}*ScaleFactor - Delta_X_stl).*SF_CAD;
% plotData{2} = (model{2}*ScaleFactor - Delta_Y_stl).*SF_CAD;
% plotData{3} = (model{3}*ScaleFactor - Delta_Z_stl).*SF_CAD;
% meshData_ac = plotData;
%
% % Colors
% C_ac(:,:,1) = color_ac(1)*ones(size(meshData_ac{1}));
% C_ac(:,:,2) = color_ac(2)*ones(size(meshData_ac{1}));
% C_ac(:,:,3) = color_ac(3)*ones(size(meshData_ac{1}));


% AC_generated

if OUTPUT_read_XLSX.Fuselage_flags.AC_original == 1   
    % Calculates the Aircraft Geometry
    model = stl2matlab(OUTPUT_read_XLSX.Fuselage_flags.STL_ac);
    ScaleFactor = OUTPUT_read_XLSX.AC_Data_flags.SF; % CAD designs are in mm, and changing to m
    SF_CAD = OUTPUT_read_XLSX.Fuselage_flags.SF_CAD;
    SF_CAD_AC = OUTPUT_read_XLSX.Fuselage_flags.SF_CAD_AC;

    if OUTPUT_read_XLSX.Fuselage_flags.STL_reading_test == 1
        max_x_ac = max(max(model{1}));
        min_x_ac = min(min(model{1}));
        max_y_ac = max(max(model{2}));
        min_y_ac = min(min(model{2}));
        max_z_ac = max(max(model{3}));
        min_z_ac = min(min(model{3}));

        Delta_X_stl = min_x_ac;
        Delta_Y_stl = -(max_y_ac - (max_y_ac - min_y_ac)/2)/2;
        Delta_Z_stl = 0;

    else
        [points_ac,triangles,tri norms] = import_stl_fast(OUTPUT_read_XLSX.Fuselage_flags.STL_ac,1);
        max_x_ac = max(points_ac(:,1));
        min_x_ac = min(points_ac(:,1));
        max_y_ac = max(points_ac(:,2));
        min_y_ac = min(points_ac(:,2));
        max_z_ac = max(points_ac(:,3));
        min_z_ac = min(points_ac(:,3));

        Delta_X_stl = points_ac(1,1);
        Delta_Y_stl = points_ac(1,2);
        Delta_Z_stl = points_ac(1,3);
    end

    % Offsets
    % plotData{1} = (model{1}*ScaleFactor + Geo_tier.x_offset_CAD - min_x_ac).*SF_CAD;
    % plotData{2} = (model{2}*ScaleFactor - (max_y_ac+min_y_ac)/2).*SF_CAD;
    % plotData{3} = (model{3}*ScaleFactor + Geo_tier.z_offset_CAD - Delta_Z_stl).*SF_CAD;
    % meshData_ac = plotData;
    
    plotData{1} = (model{1}*ScaleFactor - Delta_X_stl).*SF_CAD;
    plotData{2} = (model{2}*ScaleFactor - Delta_Y_stl).*SF_CAD;
    plotData{3} = (model{3}*ScaleFactor - Delta_Z_stl).*SF_CAD;

    plotData{1} = (model{1}*ScaleFactor - Delta_X_stl).*SF_CAD_AC;
    plotData{2} = (model{2}*ScaleFactor - Delta_Y_stl).*SF_CAD_AC;
    plotData{3} = (model{3}*ScaleFactor - Delta_Z_stl).*SF_CAD_AC;
    meshData_ac = plotData;
    
    % Colors
    C_ac(:,:,1) = color_ac(1)*ones(size(meshData_ac{1}));
    C_ac(:,:,2) = color_ac(2)*ones(size(meshData_ac{1}));
    C_ac(:,:,3) = color_ac(3)*ones(size(meshData_ac{1}));
    
    figure(Fig)
    Fig = Fig + 1;
    mesh(meshData_ac{1},meshData_ac{2},meshData_ac{3},C_ac);
    case_AC = OUTPUT_read_XLSX.AC_Data_flags.case_AC;
    st = strcat(case_AC);
    title(st,'fontsize',FS)
    title('Aircrat Original CAD')
    xlabel('y (m)')
    ylabel('z (m)')
    zlabel('z (m)')   
end
 %--------------------------------------------------


 if OUTPUT_read_XLSX.Fuselage_flags.AC_STL_Compare == 1
     figure(Fig)
     Fig = Fig + 1;
     mesh(meshData_ac{1},meshData_ac{2},meshData_ac{3},C_ac);
     title(st,'fontsize',FS)
     title('Aircrat Original and SANAID Re-Construction')
     xlabel('y (m)')
     ylabel('z (m)')
     zlabel('z (m)')
     hold on

 end
 
title(st,'fontsize',FS)
title('Aircrat Original and SANAID Re-Construction')
xlabel('y (m)')
ylabel('z (m)')
zlabel('z (m)')
hold on
st = strcat(case_AC);
title(st,'fontsize',FS)

PLOTS_Mesh_MATLAB_AC


% Plots Fuselage
% Asigns color to fuselage
model = stl2matlab(OUTPUT_read_XLSX.Fuselage_flags.STL_fus);

if OUTPUT_read_XLSX.Fuselage_flags.STL_reading_test == 1
    max_x_fus = max(max(model{1}));
    min_x_fus = min(min(model{1}));
    max_y_fus = max(max(model{2}));
    min_y_fus = min(min(model{2}));
    z_max_fus = max(max(model{3}));
    z_min_fus = min(min(model{3}));
    
    Delta_X_stl_fus = min_x_fus;
    Delta_Y_stl_fus = -(max_y_fus - (max_y_fus - min_y_fus)/2)/2;
    Delta_Z_stl_fus = 0;
else
    [points_fus,triangles,tri norms] = import_stl_fast(OUTPUT_read_XLSX.Fuselage_flags.STL_fus,1);
    max_x_fus = max(points_fus(:,1));
    min_x_fus = min(points_fus(:,1));
    max_y_fus = max(points_fus(:,2));
    min_y_fus = min(points_fus(:,2));
    z_max_fus = max(points_fus(:,3));
    z_min_fus = min(points_fus(:,3));

    Delta_X_stl_fus = points_fus(1,1);
    Delta_Y_stl_fus = points_fus(1,2);
    Delta_Y_stl_fus = 0;
    Delta_Z_stl_fus = points_fus(1,3);
end

ScaleFactor = OUTPUT_read_XLSX.AC_Data_flags.SF; % CAD designs are in mm, and changing to m
SF_CAD = OUTPUT_read_XLSX.Fuselage_flags.SF_CAD;


plotData{1} = (model{1}*ScaleFactor - Delta_X_stl_fus).*SF_CAD;
plotData{2} = (model{2}*ScaleFactor - Delta_Y_stl_fus).*SF_CAD;
% plotData{3} = (model{3}*ScaleFactor + Geo_tier.z_offset_CAD).*SF_CAD;
plotData{3} = (model{3}*ScaleFactor - Delta_Z_stl_fus).*SF_CAD;
meshData_fus = plotData;

C_fus(:,:,1) = color_fus(1)*ones(size(meshData_fus{1}));
C_fus(:,:,2) = color_fus(2)*ones(size(meshData_fus{1}));
C_fus(:,:,3) = color_fus(3)*ones(size(meshData_fus{1}));
% Plots fuselage

mesh(meshData_fus{1},meshData_fus{2},meshData_fus{3},C_fus);
case_AC = OUTPUT_read_XLSX.AC_Data_flags.case_AC;
st = strcat(case_AC);
title(st,'fontsize',FS)
title('Aircrat Original and SANAID Re-Construction')
xlabel('y (m)')
ylabel('z (m)')
zlabel('z (m)')
hold on

% Initializes
MESH_AC.dummy = 0;
% Aircraft type
switch AC_type
    case 1 % AC_type = 1 - flying wing
        %--------------------------- WING 1 ------------------------------
        MESH_AC = plot_w1(PLOTTING_UAV,color_w1,MESH_AC);
        %--------------------------- VTP ---------------------------------
        MESH_AC = plot_VTP(PLOTTING_UAV,color_vtp,MESH_AC,AC_CONFIGURATION,Geo_tier);        

    case 2 % AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
        %--------------------------- WING 1 ------------------------------
        MESH_AC = plot_w1(PLOTTING_UAV,color_w1,MESH_AC);
        %--------------------------- HTP ------------------------------
        MESH_AC = plot_HTP(PLOTTING_UAV,color_HTP,MESH_AC);
        %--------------------------- VTP ---------------------------------
        MESH_AC = plot_VTP(PLOTTING_UAV,color_vtp,MESH_AC,AC_CONFIGURATION,Geo_tier);        

    case 3 % AC_type = 3 - 3 surface: cannard + wing + HTP + VTP
        %--------------------------- WING 1 ------------------------------
        MESH_AC = plot_w1(PLOTTING_UAV,color_w1,MESH_AC);
        %--------------------------- CANARD ---------------------------------
        MESH_AC = plot_can(PLOTTING_UAV,color_HTP,MESH_AC);
        %--------------------------- HTP ------------------------------
        MESH_AC = plot_HTP(PLOTTING_UAV,color_w2,MESH_AC);
        %--------------------------- VTP ---------------------------------
        MESH_AC = plot_VTP(PLOTTING_UAV,color_vtp,MESH_AC,AC_CONFIGURATION,Geo_tier);        

    case 4 % AC_type = 4 - 2 surface: wing + V-tail
        %--------------------------- WING 1 ------------------------------
        MESH_AC = plot_w1(PLOTTING_UAV,color_w1,MESH_AC);
        %--------------------------- Vee ------------------------------
        MESH_AC = plot_vee(PLOTTING_UAV,color_vee,MESH_AC);

    case 5 % AC_type = 5 - 3 surface: cannard + wing + V-tail
        %--------------------------- WING 1 ------------------------------
        MESH_AC = plot_w1(PLOTTING_UAV,color_w1,MESH_AC);
        %--------------------------- CANARD ---------------------------------
        MESH_AC = plot_can(PLOTTING_UAV,color_can,MESH_AC);
        %--------------------------- Vee ------------------------------
        MESH_AC = plot_vee(PLOTTING_UAV,color_vee,MESH_AC);
        
    case 6 % AC_type = 6 - 2 surface: cannard + wing + VTP
        %--------------------------- WING 1 ------------------------------
        MESH_AC = plot_w1(PLOTTING_UAV,color_w1,MESH_AC);
        %--------------------------- CANARD ---------------------------------
        MESH_AC = plot_can(PLOTTING_UAV,color_can,MESH_AC);
        %--------------------------- VTP ---------------------------------
        MESH_AC = plot_VTP(PLOTTING_UAV,color_vtp,MESH_AC,AC_CONFIGURATION,Geo_tier);  


    case 7 % AC_type = 7 - 3 surface: wing + V-tail + VTP
        %--------------------------- WING 1 ------------------------------
        MESH_AC = plot_w1(PLOTTING_UAV,color_w1,MESH_AC);
        %--------------------------- Vee ------------------------------
        MESH_AC = plot_vee(PLOTTING_UAV,color_vee,MESH_AC);
        %--------------------------- VTP ---------------------------------
        MESH_AC = plot_VTP(PLOTTING_UAV,color_vtp,MESH_AC,AC_CONFIGURATION,Geo_tier);  

    case 8 % AC_type = 8 - 2 surface: wing + V-tail+V-tail

        %--------------------------- WING 1 ------------------------------
        MESH_AC = plot_w1(PLOTTING_UAV,color_w1,MESH_AC);
        %--------------------------- Vee ------------------------------
        MESH_AC = plot_vee(PLOTTING_UAV,color_vee,MESH_AC);
        %--------------------------- Vee2 ---------------------------------
        MESH_AC = plot_vee2(PLOTTING_UAV,color_vee2,MESH_AC);

end

% Engine Configuration
% Engine location
switch Engine_loc    
    case 1 % Engine_loc = 1 - under wings
        if n_eng < 2
            disp('For engines at the wing tips 2 engines need to be defined. Please do correct Input Data');
            pause
        end
       if n_eng == 2
           % Engine
           plot_2engine1(OUTPUT_read_XLSX,PLOTTING_UAV,color_eng);
           % Nacelle
           plot_2nacelle1(OUTPUT_read_XLSX,PLOTTING_UAV,color_nac);
           % Propeller disk
           plot_2prop1(OUTPUT_read_XLSX,PLOTTING_UAV,color_prop)
        elseif n_eng == 4
            % Engine and Nacelle first pair
           % Engine
           plot_2engine1(OUTPUT_read_XLSX,PLOTTING_UAV,color_eng);
           % Nacelle
           plot_2nacelle1(OUTPUT_read_XLSX,PLOTTING_UAV,color_nac);
           % Propeller disk
           plot_2prop1(OUTPUT_read_XLSX,PLOTTING_UAV,color_prop)
           % Engine and Nacelle second pair
           % Engine
           plot_2engine2(OUTPUT_read_XLSX,PLOTTING_UAV,color_eng);
           % Nacelle
           plot_2nacelle2(OUTPUT_read_XLSX,PLOTTING_UAV,color_nac);
           % Propeller disk
           plot_2prop2(OUTPUT_read_XLSX,PLOTTING_UAV,color_prop);
        end
    case 2 % Engine_loc = 2 - fuselage front

           % Engine
           plot_1engine1(OUTPUT_read_XLSX,PLOTTING_UAV,color_eng);
           % Nacelle
           plot_1nacelle1(OUTPUT_read_XLSX,PLOTTING_UAV,color_nac);
           % Propeller disk
           plot_1prop1(OUTPUT_read_XLSX,PLOTTING_UAV,color_prop)

    case 3 % Engine_loc = 3 - fuselage rear

           % Engine
           plot_1engine1(OUTPUT_read_XLSX,PLOTTING_UAV,color_eng);
           % Nacelle
           plot_1nacelle1(OUTPUT_read_XLSX,PLOTTING_UAV,color_nac);
           % Propeller disk
           plot_1prop1(OUTPUT_read_XLSX,PLOTTING_UAV,color_prop)

    case 4 % Engine_loc = 4 - wingtips
        if n_eng < 2
            disp('For engines at the wing tips 2 engines need to ve defined. Please do correct Input Data');
            pause
        end
        % Engine
        plot_2engine1(OUTPUT_read_XLSX,PLOTTING_UAV,color_eng);
        % Nacelle
        plot_2nacelle1(OUTPUT_read_XLSX,PLOTTING_UAV,color_nac);
        % Propeller disk
        plot_2prop1(OUTPUT_read_XLSX,PLOTTING_UAV,color_prop)

    case 5 % Engine_loc = 5 - wingtips 2 pairs
        if n_eng < 4
            disp('For engines at the wing tips 2 engines need to be defined. Please do correct Input Data');
            pause
        end
        if n_eng == 2
            % Engine
            plot_2engine1(OUTPUT_read_XLSX,PLOTTING_UAV,color_eng);
            % Nacelle
            plot_2nacelle1(OUTPUT_read_XLSX,PLOTTING_UAV,color_nac);
            % Propeller disk
            plot_2prop1(OUTPUT_read_XLSX,PLOTTING_UAV,color_prop)
        elseif n_eng == 4
            % Engine and Nacelle first pair
            % Engine
            plot_2engine1(OUTPUT_read_XLSX,PLOTTING_UAV,color_eng);
            % Nacelle
            plot_2nacelle1(OUTPUT_read_XLSX,PLOTTING_UAV,color_nac);
            % Propeller disk
            plot_2prop1(OUTPUT_read_XLSX,PLOTTING_UAV,color_prop)
            % Engine and Nacelle second pair
            % Engine
            plot_2engine2(OUTPUT_read_XLSX,PLOTTING_UAV,color_eng);
            % Nacelle
            plot_2nacelle2(OUTPUT_read_XLSX,PLOTTING_UAV,color_nac);
            % Propeller disk
            plot_2prop2(OUTPUT_read_XLSX,PLOTTING_UAV,color_prop);
        end
end
grid on

if OUTPUT_read_XLSX.Fuselage_flags.AC_STL_Compare == 1
    
    % Calculates the Aircraft Geometry
    model = stl2matlab(OUTPUT_read_XLSX.Fuselage_flags.STL_ac);
    ScaleFactor = OUTPUT_read_XLSX.AC_Data_flags.SF; % CAD designs are in mm, and changing to m
    SF_CAD = OUTPUT_read_XLSX.Fuselage_flags.SF_CAD;
    SF_CAD_AC = OUTPUT_read_XLSX.Fuselage_flags.SF_CAD_AC;
    
    if OUTPUT_read_XLSX.Fuselage_flags.STL_reading_test == 1
        max_x_ac = max(max(model{1}));
        min_x_ac = min(min(model{1}));
        max_y_ac = max(max(model{2}));
        min_y_ac = min(min(model{2}));
        max_z_ac = max(max(model{3}));
        min_z_ac = min(min(model{3}));

        Delta_X_stl = min_x_ac;
        Delta_Y_stl = -(max_y_ac - (max_y_ac - min_y_ac)/2)/2;
        Delta_Z_stl = 0;

    else
        [points_ac,triangles,tri norms] = import_stl_fast(OUTPUT_read_XLSX.Fuselage_flags.STL_ac,1);
        max_x_ac = max(points_ac(:,1));
        min_x_ac = min(points_ac(:,1));
        max_y_ac = max(points_ac(:,2));
        min_y_ac = min(points_ac(:,2));
        max_z_ac = max(points_ac(:,3));
        min_z_ac = min(points_ac(:,3));

        Delta_X_stl = points_ac(1,1);
        Delta_Y_stl = points_ac(1,2);
        Delta_Z_stl = points_ac(1,3);
    end
     
    
    % 
    % Offsets
    % plotData{1} = (model{1}*ScaleFactor + Geo_tier.x_offset_CAD - min_x_ac).*SF_CAD;
    % plotData{2} = (model{2}*ScaleFactor - (max_y_ac+min_y_ac)/2).*SF_CAD;
    % plotData{3} = (model{3}*ScaleFactor + Geo_tier.z_offset_CAD - Delta_Z_stl).*SF_CAD;
    % meshData_ac = plotData;
    
    plotData{1} = (model{1}*ScaleFactor - Delta_X_stl).*SF_CAD_AC;
    plotData{2} = (model{2}*ScaleFactor - Delta_Y_stl).*SF_CAD_AC;
    plotData{3} = (model{3}*ScaleFactor - Delta_Z_stl).*SF_CAD_AC;
    meshData_ac = plotData;
    
    % Colors
    C_ac(:,:,1) = color_ac(1)*ones(size(meshData_ac{1}));
    C_ac(:,:,2) = color_ac(2)*ones(size(meshData_ac{1}));
    C_ac(:,:,3) = color_ac(3)*ones(size(meshData_ac{1}));
    
    mesh(meshData_ac{1},meshData_ac{2},meshData_ac{3},C_ac);
end
hold off


%% Defines limits of plots axis
z_max_fus = PLOTTING_UAV.z_max_fus;
z_min_fus = PLOTTING_UAV.z_min_fus;
% Creates Margin to all sizes of plots around the airplane
Delta_Plot = 0.05;
b_w1 = Geo_tier.b_w1;
MAX_X = L_fus + Delta_Plot;
MAX_Y = b_w1 + Delta_Plot;
MAX_Z = z_max_fus + Delta_Plot;
x_min_PLOT = -Delta_Plot;
x_max_PLOT = MAX_X;
y_min_PLOT = -MAX_Y/2 ;
y_max_PLOT = MAX_Y/2 ;
z_min_PLOT = (z_min_fus);
z_max_PLOT = MAX_Z;
axis([x_min_PLOT x_max_PLOT y_min_PLOT y_max_PLOT z_min_PLOT z_max_PLOT]);
axis equal

PLOTTING_DATA = PLOTTING_UAV;
if SAVE_FIGS == 1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('Geometry');
    name   = strcat(prefix,st);
    % fname = OUTPUT_read_XLSX.PLOT_flags.fname;
            % Unified plot Options             
            SAVE_types(fname,name,gca,gcf); 
%             saveas(gcf,fullfile(fname, name),'fig');
%             saveas(gcf,fullfile(fname, name),'pdf');
%             saveas(gcf,fullfile(fname, name),'bmp');
%             saveas(gcf,fullfile(fname, name),'png');

end

if Video_3D == 1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('VID_3D_');
    name1   = strcat(prefix,st);
    VID_TXT = name1;
    OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
    CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10],VID_TXT,OptionZ)
end
% Fig = Fig + 1;