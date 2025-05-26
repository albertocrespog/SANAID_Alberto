close all
clear all
% Defines Generic Path so that files can be organiez in folders
% change for the different versions
addpath(genpath('../MATLAB'))
% addpath(genpath('../ProyectoEMERGENTIA/XFLR5'))
% addpath(genpath('../EMERGENTIA/XFLR5'))
% Units conversion
conv_UNITS = conversion_UNITS;

R2D = conv_UNITS.R2D;
D2R = conv_UNITS.D2R;

% Initial flight conditions
h = 0;
% Cruise speed
V = 25;
V_max_CR = 1.25*V; % Max Speed 

% Flight Conditons
[Temp,rho,p,a]=atmos_inter(h);
Mach = V/a;

% Saves data regarding flight conditions
V_Performance.h = h;
V_Performance.Mach = Mach;
V_Performance.Temp = Temp;
V_Performance.rho = rho;
V_Performance.p = p;
V_Performance.a = a;
V_Performance.V_CR = V;
V_Performance.V_max_CR = V_max_CR;

Data_ATM.Temp = Temp;
Data_ATM.rho = rho;
Data_ATM.p = p;
Data_ATM.a = a;

% Initial Geometry
Geo_tier = Geometric_Data_tier;

m_T0 = 18.8;
Weights.m_T0 = m_T0;

%% Selection of TXT that are used for the aerodynamic analysis
% Select txt files associated with each aerodynamic study
% Make sure to check with case asssigned in read_aero_files_Aug2018.m
index_w1 = 3; % Front Wing W1
index_w2 = 9; % Rear Wing W2

% Allows the user to select the min AoA used to determine Curve Lift Slope
alpha_selected_w1 = 0; % (degs)
alpha_selected_w2 = 0; % (degs)
Design_criteria.index_w1 = index_w1;
Design_criteria.index_w2 = index_w2;
Design_criteria.alpha_selected_w1 = alpha_selected_w1;
Design_criteria.alpha_selected_w2 = alpha_selected_w2;

% Design choices for the incidence of the different surfaces
Design_criteria.i_w1 = 4*D2R; % incidence of Front Wing 
Design_criteria.i_w2 = 0*D2R; % incidence of Rear Wing

% Flight Safe Margin to calculate aerodynamic properties
Flight_SF = 1.2; % Stall Safe Margin Conditions
Design_criteria.Flight_SF = Flight_SF;

%% Trimmed conditions
% Obtains Trimmed conditions and Derivatives
Data_Trim = Trim_Cefiro_1214(h,V);
% Cefiro Derivatives need to be changed
Data_Der = get_derivatives_Cefiro_1214(Data_Trim,h,V);

% Routine that generates States
[x,delta]=Generate_States(Data_Trim,Data_Der,h,V);

D = Geo_tier.D_prop;
alpha = Data_Trim.alpha_1;
% Number of Engine used
k = 4;
% Initializes figures
fig = 0;

%% Obtains propulsion values
% Desired Force (Newtons)
Fdes = 30;
% Gets the engine properties (for each engines)
[Propulsion] = get_EngineProperties(Data_Trim,V,alpha,Geo_tier,k,Fdes);

%%
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
l_fus2 = 1.6797;
Ratio_SCALE = l_fus2/l_fus1;

% Being used only if STL model used to determine geometry of fuselage
% Nth = size(plotData{1});
% STLnSections = Nth(1);% number of sectins (x coordinate) (eliminated 1 for convergence)
% STLnPoints = Nth(2); % number of elements per section (y and z coordinates)

% Determines fuselage information
[Body_Geo,meshData] = Calc_Body_Geometry_Dec2018(XFLR5_file,STL_file,Geo_tier, nPoints, nSections,conversion_UNITS,lecture_Geo,Ratio_SCALE);

% Genera fuselage from STL of CAD
STL_PLOT = 0;
if STL_PLOT == 1
    % Loads STL file into a mesh
    model = stl2matlab('ProVant.stl');
    ScaleFactor = 1/1000; % CAD designs are in mm, and changing to m
    plotData{1} = model{1}*ScaleFactor;
    plotData{2} = model{2}*ScaleFactor;
    plotData{3} = model{3}*ScaleFactor;
    % Saves the mesh data from the STL CAD graphics
    meshData = plotData;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SAVE_FIGS = 0;  % Save figs with different formats
PRINT_PLOTS_XFLR5 = 0; % Prints plots - Aero
PRINT_PLOTS_XFLR5_POLAR = 1; % Prints plots - Aero
PRINT_PLOTS_AERO = 0; % Prints plots aerodynamic analysis
PRINT_PLOTS_3D = 0; % Prints plots for 3D
PRINT_PLOTS_XAC = 0; % print Plots for Stimation ofXAC
PRINT_PLOTS_CT_Model = 0; % print Plots for Propulsive Models
figures_CT_model = 0; % prints plots of propulsive model
% size of plot letter
LS = 1; % Line size
FS = 12; % Text Font size
LFS = 8; % Legend Font Size
Fig = 0;
Video_3D = 0; % Saves video 3D

% Determines the plots that want to show
% VECTOR = [1,2,3]; % Compares Aerodynamic properties with moments meassured on Leading Edge and in the Aerodynamic Center for W_1

% The VECTOR in order to analyze the results of XFLR consist of 4 grous of
% vector each for: w1, w2, vtp, and full airplane
VECTOR_XFLR5.v{1} = [1,2,3,4,5,6,7,8]; % compares VLM &/or LLT for w1
% VECTOR_XFLR5.v{2} = [9,10,11,12,13,14]; % compares VLM & LLT for w2
VECTOR_XFLR5.v{2} = [9,10]; % compares VLM & LLT for w2

% Desices which cases are compared
% compare = 1 - w1 wing alone
% compare = 2 - w2 wing alone
% compare = 3 - w1 & w2 
compare = 2;
VECTOR_XFLR5.compare = compare; % compares VLM & LLT for w2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stores information for the plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data for the legend
% V_CR = V_Performance.V_CR; % Cruise speed
% h_i = V_Performance.h_i; % Altitude
% Generates the Legend
[mark_Type] = Generates_plot_Info(V,h); % Generates the style of the lines and the Legend
% Read aerodnamic data
[casos prefix mark_legend X_OC] = read_aero_files_Aug2018(V_Performance);
[DATA_Ae] = Read_data_Aero(casos); 

Plot_Options.LS = LS;
Plot_Options.FS = FS;
Plot_Options.LFS = LFS;
Plot_Options.Fig = Fig;
Plot_Options.Video_3D = Video_3D;
Plot_Options.SAVE_FIGS = SAVE_FIGS;
Plot_Options.PRINT_PLOTS_XFLR5 = PRINT_PLOTS_XFLR5;
Plot_Options.PRINT_PLOTS_AERO = PRINT_PLOTS_AERO;
Plot_Options.PRINT_PLOTS_3D = PRINT_PLOTS_3D;
Plot_Options.PRINT_PLOTS_XFLR5_POLAR = PRINT_PLOTS_XFLR5_POLAR;
Plot_Options.PRINT_PLOTS_XAC = PRINT_PLOTS_XAC;
Plot_Options.PRINT_PLOTS_CT_Model = PRINT_PLOTS_CT_Model;
Plot_Options.mark_Type = mark_Type;
Plot_Options.mark_legend = mark_legend;
% Plot_Options.VECTOR = VECTOR;
Plot_Options.VECTOR_XFLR5 = VECTOR_XFLR5;
Plot_Options.figures_CT_model = figures_CT_model;
PL = 1;

%% Theoretical Aerodynamics
% Fusion Aerodynamic Properties
[Aero_TH] = Polar_Literatura(h,V,Geo_tier,conv_UNITS,Body_Geo);

%% Theoretical Aerodynamics
% Define how many cases have been analyzed to determine Xac
Casos_XAC_w1 = 6;
Casos_XAC_w2 = 6;

% Defines Degres that are used to determine XAC
degrees_XAC = [2.5 5 7.5 10 12.5 15 17.5];
XAC_vec = [0 0.05 0.10 0.15 0.20 0.25];

CASOS.Casos_XAC_w1 = Casos_XAC_w1;
CASOS.Casos_XAC_w2 = Casos_XAC_w2;
% mean aerodynamic chord according to XFLR5
MAC_XFLR5 = 0.205;
S_w1_XFLR5 = 0.452;
XNP_XFLR5 = -0.031;
XFLR5_DATA.MAC_XFLR5 = MAC_XFLR5;
XFLR5_DATA.S_w1_XFLR5 = S_w1_XFLR5;
XFLR5_DATA.XNP_XFLR5 = XNP_XFLR5;

% Etimation of the MAC
MAC_Stimation = 0;
% Premiliminary Aerodynamic Design
[Aero,DATA_PL,Fig] = Aerodynamic_Design_Aug2018(Geo_tier,...
    Weights,conv_UNITS,Design_criteria,V_Performance,CASOS,degrees_XAC,fig,XFLR5_DATA,...
    MAC_Stimation,DATA_Ae,X_OC,Plot_Options,VECTOR_XFLR5);

%% 
% Generation of Forces 
% Wind Components
Wx = 0;
Wy = 0;
Wz = 0;

Wind.Wx = Wx;
Wind.Wy = Wy;
Wind.Wz = Wz;

[F, M,DynVar] = Get_ForcesMoments(x,delta,Data_Trim,Data_Der,h,V,Wind,Data_ATM,Propulsion)

n = Propulsion.n;
%% 
% Plots graphics
if PRINT_PLOTS_XFLR5 == 1
    [Fig] = plot_XFLR5_Aero_Aug2018(Aero,DATA_Ae,DATA_PL,Plot_Options,Fig,Design_criteria,prefix);
elseif PRINT_PLOTS_XFLR5_POLAR == 1
    [Fig] = plot_XFLR5_Polar_Aug2018(Aero,DATA_Ae,DATA_PL,Plot_Options,Fig,Design_criteria,prefix)
elseif PRINT_PLOTS_AERO == 1
    [Fig] = plot_Aero_Aug2018(Aero,DATA_Ae,DATA_PL,V_Performance,Plot_Options,Fig,prefix);
elseif PRINT_PLOTS_XAC == 1
    [XAC,Fig] = plot_Generate_XAC(DATA_Ae,Aero,DATA_PL,CASOS,degrees_XAC,Fig,...
    XFLR5_DATA,V_Performance,X_OC,Plot_Options,VECTOR_XFLR5)
elseif PRINT_PLOTS_CT_Model == 1 % Generates Propulsion Plots
    [Fig] = Generates_Plots_PropulsionModels(Propulsion,Data_Trim,V,alpha,n,D,k,Plot_Options,Fig);
elseif PRINT_PLOTS_3D == 1 % Generates Data for Ploting 3D Mesh
    [Fig] = plot_GEOMETRY_2018(Geo_tier,Plot_Options,Body_Geo,meshData,Propulsion,Fig);
end

