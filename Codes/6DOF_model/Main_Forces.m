close all
clear all
% Defines Generic Path so that files can be organized in folders
% change for the different versions
addpath(genpath('../MATLAB'))

% Units conversion
conv_UNITS = conversion_UNITS;

% Initial flight conditions
h = 200; % altitude [m]
% Cruise speed
V = 25; % [m/s]
V_max_CR = 1.25*V; % Max Speed [m/s]

% Flight Conditons obtained using Standard Atmospheric Properties
% [Temp,rho,p,a]=atmos_inter(h)
Toffset = 0;
[rho,a,Temp,p,kvisc,ZorH]=stdatmo(h,Toffset,'SI');
Mach = V/a; % Mach number

% Saves data regarding flight conditions
V_Performance.h = h;
V_Performance.Mach = Mach;
V_Performance.Temp = Temp;
V_Performance.rho = rho;
V_Performance.p = p;
V_Performance.a = a;
V_Performance.V_CR = V;
V_Performance.V_max_CR = V_max_CR;

% Atmospheric properties
%           rho:   Density            kg/m^3          slug/ft^3
%           a:     Speed of sound     m/s             ft/s
%           T:     Temperature        K              R
%           P:     Pressure           Pa              lbf/ft^2
%           nu:    Kinem. viscosity   m^2/s           ft^2/s
%           ZorH:  Height or altitude m               ft
Data_ATM.rho = rho;
Data_ATM.a = a;
Data_ATM.Temp = Temp;
Data_ATM.p = p;
Data_ATM.kvisc = kvisc;
Data_ATM.ZorH = ZorH;

% Initial Geometry
Geo_tier = Geometric_Data_tier;




% Obtains Trimmed conditions and Derivatives
Data_Trim = Trim_Cefiro_1214(h,V);
% Cefiro Derivatives need to be changed
Data_Der = get_derivatives_Cefiro_1214(Data_Trim,h,V);

% Routine that generates States
[x,delta]=Generate_States(Data_Trim,Data_Der,h,V);

D = Geo_tier.D_prop;
% n = 5000/60;
alpha = Data_Trim.alpha_1;
% Number of Prop used
k = 4;
% Initializes figures
fig = 0;

% Desired Force (Newtons)
Fdes = 30;
% Gets the engine properties (for each engines)
[Propulsion] = get_EngineProperties(Data_Trim,V,alpha,Geo_tier,k,Fdes);

% Cruise speed
V_CR = 25;
V_max_CR = 1.25*V_CR; % Max Speed 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates Geometry Geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XFLR5_file = strcat('fuse_pepi_desplazado.txt');
XFLR5_file = strcat('Fuselage_DATA_scaled.txt');
STL_file = strcat('Fus_ProVant_v4.0.stl');

% Pepiño
%--------------------------- FUSELAJE ---------------------------------
nSections = 96;% number of sectins (x coordinate) (eliminated 1 for convergence)
nPoints = 67; % number of elements per section (y and z coordinates)
lecture_Geo = 5; % lineas a partir desde donde empieza a leer
% EMERGENTIA
%--------------------------- FUSELAJE ---------------------------------
nSections = 9;% number of sectins (x coordinate) (eliminated 1 for convergence)
nPoints = 50; % number of elements per section (y and z coordinates)
lecture_Geo = 5; % lineas a partir desde donde empieza a leer
%--------------------------- CALCULATES GEOMETRY DATA ---------------------------------
% Scale from fuselage in XFLR5 to CAD
l_fus1 = 0.8799;
l_fus2 = Geo_tier.l_fus;
Ratio_SCALE = l_fus2/l_fus1;

% Being used only if STL model used to determine geometry of fuselage
% Nth = size(plotData{1});
% STLnSections = Nth(1);% number of sectins (x coordinate) (eliminated 1 for convergence)
% STLnPoints = Nth(2); % number of elements per section (y and z coordinates)

% Determines fuselage information
[Body_Geo,meshData] = Calc_Body_Geometry_Dec2018(XFLR5_file,STL_file,Geo_tier, nPoints, nSections,conversion_UNITS,lecture_Geo,Ratio_SCALE);

% Genera fuselage from STL of CAD
STL_PLOT = 1;
if STL_PLOT == 1
    % Loads STL file into a mesh
    model = stl2matlab('ProVant.stl');
%     model = stl2matlab('VTOL_Single.stl');
    ScaleFactor = 1/1000; % CAD designs are in mm, and changing to m
    plotData{1} = model{1}*ScaleFactor + Geo_tier.x_offset_CAD;
    plotData{2} = model{2}*ScaleFactor;
    plotData{3} = model{3}*ScaleFactor + Geo_tier.z_offset_CAD;
    % Saves the mesh data from the STL CAD graphics
    meshData = plotData;
end

% Fusion Aerodynamic Properties
[Aero_TH] = Polar_Literatura(h,V,Geo_tier,conv_UNITS,Body_Geo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SAVE_FIGS = 0;  % Save figs with different formats
PRINT_PLOTS_XFLR5 = 0; % Prints plots
PRINT_PLOTS_AERO = 0; % Prints plots
PRINT_PLOTS_PROPULSION = 0 ; % Prints plots
figures_CT_model = 1; % prints plots of propulsive model
% size of plot letter
LS = 2; % Line size
FS = 12; % Text Font size
LFS = 8; % Legend Font Size
Fig = 0;
Video_3D = 0; % Saves video 3D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stores information for the plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data for the legend
% V_CR = V_Performance.V_CR; % Cruise speed
% h_i = V_Performance.h_i; % Altitude
% Generates the Legend
[mark_Type] = Generates_plot_Info(V,h); % Generates the style of the lines and the Legend

% Determines the plots that want to show
% VECTOR = [1,2,5,6,9,10]; % Defines the number of Aerodynamic Cases to be analyzed: See "read_aero_files_Aug2018.m" to select them
VECTOR = [2,5]; % Compares Aerodynamic properties with moments meassured on Leading Edge and in the Aerodynamic Center for W_1
% VECTOR = [5,7]; % Compares Aerodynamic properties with moments meassured on Leading Edge and in the Aerodynamic Center for W_2
% VECTOR = [9,10]; % Compares Aerodynamic properties with moments meassured on Leading Edge and in the Aerodynamic Center for VTP
% VECTOR = [6,7,8]; % Compares Aerodynamic properties with moments meassured on Leading Edge and in the Aerodynamic Center for W_2

% The VECTOR in order to analyze the results of XFLR consist of 4 grous of
% vector each for: w1, w2, vtp, and full airplane
VECTOR_XFLR5.v{1} = [3,5]; % compares VLM & LLT for w1
VECTOR_XFLR5.v{2} = [7,10]; % compares VLM & LLT for w2
VECTOR_XFLR5.v{3} = [12,14]; % compares VLM & LLT for vtp
VECTOR_XFLR5.v{4} = [15]; % compares VLM & LLT for full airplane

Plot_Options.LS = LS;
Plot_Options.FS = FS;
Plot_Options.LFS = LFS;
Plot_Options.Fig = Fig;
Plot_Options.Video_3D = Video_3D;
Plot_Options.SAVE_FIGS = SAVE_FIGS;
Plot_Options.PRINT_PLOTS_XFLR5 = PRINT_PLOTS_XFLR5;
Plot_Options.PRINT_PLOTS_AERO = PRINT_PLOTS_AERO;
Plot_Options.PRINT_PLOTS_PROPULSION = PRINT_PLOTS_PROPULSION;
Plot_Options.mark_Type = mark_Type;
% Plot_Options.mark_legend = mark_legend;
Plot_Options.VECTOR = VECTOR;
Plot_Options.VECTOR_XFLR5 = VECTOR_XFLR5;
Plot_Options.figures_CT_model = figures_CT_model;
PL = 1;

n = Propulsion.n;

% Generates Propulsion Plots
[Fig] = Generates_Plots_PropulsionModels(Propulsion,Data_Trim,V,alpha,n,D,k,Plot_Options,Fig);

% Generates Data for Ploting 3D Mesh
[Fig] = plot_GEOMETRY_2018(Geo_tier,Plot_Options,Body_Geo,meshData,Propulsion,Fig);

% Wind Components
Wx = 0;
Wy = 0;
Wz = 0;

Wind.Wx = Wx;
Wind.Wy = Wy;
Wind.Wz = Wz;

[F, M,DynVar] = Get_ForcesMoments(x,delta,Data_Trim,Data_Der,h,V,Wind,Data_ATM,Propulsion,Geo_tier);

