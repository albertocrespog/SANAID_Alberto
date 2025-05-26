function Fig = plot_GEOMETRY_2022_v2(Geo_tier,Plot_Options,Body_Geo,meshData,Prop_data,Fig,AC_CONFIGURATION,COLOR_scheme,case_AC,Performance,OUTPUT_read_XLSX,filenameS)

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

% Stores the names of folders zand address
% filename = 'Results\XX_YYYYY\';
% filenamed = 'DATA';
% filenamec = 'Figs';
% filename_DATA = strcat(filename,filenamed);
% filename_Plots = strcat(filename,filenamec);
fname = filenameS.plots;

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

if OUTPUT_read_XLSX.Fuselage_flags.Use_Storing_DATA == 0
    %% Generation of Geometry
    [Fig,Body_Geo,plotData] = Generate_Plots_Body_Geometry(XFLR5_file,STL_files,Geo_tier,nPoints, nSections,conversion_UNITS,lecture_Geo,Ratio_SCALE,SF,OUTPUT_read_XLSX,Fig,filenameS);
    Saving_data_PLOT_FUS(Body_Geo,plotData,OUTPUT_read_XLSX,filenameS);
else
    [Body_Geo,plotData] = Read_data_PLOT_FUS(OUTPUT_read_XLSX,filenameS);
    Fig = 1;
end

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

% AC_generated
if OUTPUT_read_XLSX.Fuselage_flags.AC_original == 1  
    figure(Fig)
    Fig = Fig + 1;
    PLOTS_Real_AC(OUTPUT_read_XLSX,PLOTTING_UAV,COLOR_scheme,Body_Geo);
    case_AC = OUTPUT_read_XLSX.AC_Data_flags.case_AC;
    st = strcat(case_AC);
    title(st,'fontsize',FS)
    title('Aircrat Original CAD')
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
    if SAVE_FIGS == 1
        prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
        st = strcat('Geometry_CAD');
        name   = strcat(prefix,st);
        SAVE_types(fname,name,gca,gcf);
    end

    %% Defines limits of plots axis
    z_max_fus = PLOTTING_UAV.z_max_fus;
    z_min_fus = PLOTTING_UAV.z_min_fus;
    % Creates Margin to all sizes of plots around the airplane
    Ddelta_max = 1.25;
    Ddelta_min = 0.25;
    b_w1 = Geo_tier.b_w1;
    MAX_X = L_fus*Ddelta_max;
    MAX_Y = b_w1*Ddelta_max;
    MAX_Z = z_max_fus*Ddelta_max;
    x_min_PLOT = -L_fus*Ddelta_min;
    x_max_PLOT = MAX_X;
    y_min_PLOT = -MAX_Y/2 ;
    y_max_PLOT = MAX_Y/2 ;
    z_min_PLOT = -MAX_Z;
    z_max_PLOT = MAX_Z;
    axis([x_min_PLOT x_max_PLOT y_min_PLOT y_max_PLOT z_min_PLOT z_max_PLOT]);
    axis equal
end

if OUTPUT_read_XLSX.Fuselage_flags.AC_STL_Compare == 1
    figure(Fig)
    Fig = Fig + 1;
    % Plots AC Original
    PLOTS_Real_AC(OUTPUT_read_XLSX,PLOTTING_UAV,COLOR_scheme);
    st = strcat(case_AC);
    title(st,'fontsize',FS)
    title('Aircrat Original CAD and SANAID Re-Construction')
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
    hold on

    %   %% Defines limits of plots axis
    % z_max_fus = PLOTTING_UAV.z_max_fus;
    % z_min_fus = PLOTTING_UAV.z_min_fus;
    % % Creates Margin to all sizes of plots around the airplane
    % Ddelta_max = 1.25;
    % Ddelta_min = 0.25;
    % b_w1 = Geo_tier.b_w1;
    % MAX_X = L_fus*Ddelta_max;
    % MAX_Y = b_w1*Ddelta_max;
    % MAX_Z = z_max_fus*Ddelta_max;
    % x_min_PLOT = -L_fus*Ddelta_min;
    % x_max_PLOT = MAX_X;
    % y_min_PLOT = -MAX_Y/2 ;
    % y_max_PLOT = MAX_Y/2 ;
    % z_min_PLOT = -MAX_Z;
    % z_max_PLOT = MAX_Z;
    % axis([x_min_PLOT x_max_PLOT y_min_PLOT y_max_PLOT z_min_PLOT z_max_PLOT]);
    % axis equal

    % Plots AC Mesh MATLAB
    PLOTS_Mesh_MATLAB_AC(OUTPUT_read_XLSX,PLOTTING_UAV,COLOR_scheme,AC_type,Engine_loc,AC_CONFIGURATION,Geo_tier);
    PLOTS_Mesh_MATLAB_Propulsion(OUTPUT_read_XLSX,PLOTTING_UAV,COLOR_scheme,AC_type,Engine_loc,n_eng)
    
    hold off
    grid on

    if SAVE_FIGS == 1
        prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
        st = strcat('Geometry_CAD_SANAID');
        name   = strcat(prefix,st);
        SAVE_types(fname,name,gca,gcf);
    end

    %% Defines limits of plots axis
    z_max_fus = PLOTTING_UAV.z_max_fus;
    z_min_fus = PLOTTING_UAV.z_min_fus;
    % Creates Margin to all sizes of plots around the airplane
    Ddelta_max = 1.25;
    Ddelta_min = 0.25;
    b_w1 = Geo_tier.b_w1;
    MAX_X = L_fus*Ddelta_max;
    MAX_Y = b_w1*Ddelta_max;
    MAX_Z = z_max_fus*Ddelta_max;
    x_min_PLOT = -L_fus*Ddelta_min;
    x_max_PLOT = MAX_X;
    y_min_PLOT = -MAX_Y/2 ;
    y_max_PLOT = MAX_Y/2 ;
    z_min_PLOT = -MAX_Z;
    z_max_PLOT = MAX_Z;
    
end

if OUTPUT_read_XLSX.Fuselage_flags.AC_generated == 1
    figure(Fig)
    Fig = Fig + 1;
    st = strcat(case_AC);
    title(st,'fontsize',FS)
    title('Aircrat Generated Mesh - SANAID Re-Construction')
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
    hold on
    
    % Plots AC Mesh MATLAB
    PLOTS_Mesh_MATLAB_AC(OUTPUT_read_XLSX,PLOTTING_UAV,COLOR_scheme,AC_type,Engine_loc,AC_CONFIGURATION,Geo_tier);
    PLOTS_Mesh_MATLAB_Propulsion(OUTPUT_read_XLSX,PLOTTING_UAV,COLOR_scheme,AC_type,Engine_loc,n_eng)
    
    hold off
    grid on
    
    if SAVE_FIGS == 1
        prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
        st = strcat('Geometry_SANAID');
        name   = strcat(prefix,st);
        SAVE_types(fname,name,gca,gcf);
    end

    %% Defines limits of plots axis
    z_max_fus = PLOTTING_UAV.z_max_fus;
    z_min_fus = PLOTTING_UAV.z_min_fus;
    % Creates Margin to all sizes of plots around the airplane
    Ddelta_max = 1.25;
    Ddelta_min = 0.25;
    b_w1 = Geo_tier.b_w1;
    MAX_X = L_fus*Ddelta_max;
    MAX_Y = b_w1*Ddelta_max;
    MAX_Z = z_max_fus*Ddelta_max;
    x_min_PLOT = -L_fus*Ddelta_min;
    x_max_PLOT = MAX_X;
    y_min_PLOT = -MAX_Y/2 ;
    y_max_PLOT = MAX_Y/2 ;
    z_min_PLOT = -MAX_Z;
    z_max_PLOT = MAX_Z;
    ax = gca;
    ax.Interactions = [rotateInteraction dataTipInteraction];

end



if Video_3D == 1
    prefix = strcat(OUTPUT_read_XLSX.PLOT_flags.prefix);
    st = strcat('VID_3D_');
    name1   = strcat(prefix,st);
    VID_TXT = name1;
    OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
    CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10],VID_TXT,OptionZ)
end