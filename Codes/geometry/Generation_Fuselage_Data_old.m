function [Fig,Body_Geo,meshData] = Generation_Fuselage_Data_old(Geo_tier,OUTPUT_read_XLSX,Fig)% Defines Propulsion DATA

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
[Fig,Body_Geo,meshData] = Calc_Body_Geometry_Dec2018_v1(XFLR5_file,STL_files,Geo_tier, nPoints, nSections,conversion_UNITS,lecture_Geo,Ratio_SCALE,SF,OUTPUT_read_XLSX,Fig);
% 
% Calc_Body_Geometry_Feb2021(XFLR5_file,STL_file,Geo_tier,nPoints, nSections,conversion_UNITS,lecture_Geo,Ratio_SCALE,SF,plotData,DATA_FUS)
