close all
clear all
% Defines Generic Path so that files can be organiez in folders
% change for the different versions
addpath(genpath('../EMERGENTIA'))
% Units conversion
conv_UNITS = conversion_UNITS;

R2D = conv_UNITS.R2D;
D2R = conv_UNITS.D2R;

% Flight Conditons
h = 200;
[Temp,rho,p,a]=atmos_inter(h);
V = 25;
Mach = V/a;
% Cruise speed
V_max_CR = 1.25*V; % Max Speed 

% Saves data regarding flight conditions
V_Performance.h = h;
V_Performance.Mach = Mach;
V_Performance.Temp = Temp;
V_Performance.rho = rho;
V_Performance.p = p;
V_Performance.a = a;
V_Performance.V_CR = V;
V_Performance.V_max_CR = V_max_CR;

% Data_Trim = Trim_Cefiro_1214(h,V);
% Data_Der = get_derivatives_Cefiro_1214(Data_Trim,h,V);
% 
% % Routine that generates States
% [x,delta]=Generate_States(Data_Trim,Data_Der,h,V)
% D = 0.6043;
% n = 5000/60;
% alpha = Data_Trim.alpha_1;
% k = 4;
% % Gets the engine properties
% [Propulsion] = get_EngineProperties(Data_Trim,V,alpha,n,D,k) 

% Initial Geometry
fig = 0;
Geo_tier = Geometric_Data_tier;

m_T0 = 18.8;
Weights.m_T0 = m_T0;

% Select txt files associated with each aerodynamic study
% Make sure to check with case asssigned in read_aero_files_Aug2018.m
index_w1 = 3; % Front Wing W1
index_w2 = 6; % Rear Wing W2

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

% Define how many cases have been analyzed to determine Xac
Casos_XAC_w1 = 6;
Casos_XAC_w2 = 6;

% Defines Degres that are used to determine XAC
degrees_XAC = [-5 -2.5 0 2.5 5 7.5 10];
XAC_vec = [0 0.05 0.10 0.15 0.20 0.25];

CASOS.Casos_XAC_w1 = Casos_XAC_w1;
CASOS.Casos_XAC_w2 = Casos_XAC_w2;
% Premiliminary Aerodynamic Design
[Aero,DATA_Ae,DATA_PL,fig] = Aerodynamic_Design_Aug2018(Geo_tier,Weights,conv_UNITS,Design_criteria,V_Performance,CASOS,degrees_XAC,fig,XAC_vec)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stores information for the plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data for the legend
V_CR = V_Performance.V_CR; % Cruise speed
% Generates the Legend
[mark_Type,mark_legend] = Generates_plot_Info(V_CR,h); % Generates the style of the lines and the Legend

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SAVE_FIGS = 0;  % Save figs with different formats
PRINT_PLOTS_XFLR5 = 0; % Prints plots
PRINT_PLOTS_AERO = 0; % Prints plots
% size of plot letter
LS = 2; % Line size
FS = 12; % Text Font size
LFS = 8; % Legend Font Size
Fig = 0;
Video_3D = 0; % Saves video 3D

% Determines the plots that want to show
% VECTOR = [1,2,5,6,9,10]; % Defines the number of Aerodynamic Cases to be analyzed: See "read_aero_files_Aug2018.m" to select them
VECTOR = [1,2,3]; % Compares Aerodynamic properties with moments meassured on Leading Edge and in the Aerodynamic Center for W_1
% VECTOR = [5,7]; % Compares Aerodynamic properties with moments meassured on Leading Edge and in the Aerodynamic Center for W_2
% VECTOR = [9,10]; % Compares Aerodynamic properties with moments meassured on Leading Edge and in the Aerodynamic Center for VTP
% VECTOR = [6,7,8]; % Compares Aerodynamic properties with moments meassured on Leading Edge and in the Aerodynamic Center for W_2

% The VECTOR in order to analyze the results of XFLR consist of 4 grous of
% vector each for: w1, w2, vtp, and full airplane
VECTOR_XFLR5.v{1} = [1,2,3,4,5]; % compares VLM & LLT for w1
VECTOR_XFLR5.v{2} = [7,10]; % compares VLM & LLT for w2

Plot_Options.LS = LS;
Plot_Options.FS = FS;
Plot_Options.LFS = LFS;
Plot_Options.Fig = Fig;
Plot_Options.Video_3D = Video_3D;
Plot_Options.SAVE_FIGS = SAVE_FIGS;
Plot_Options.PRINT_PLOTS_XFLR5 = PRINT_PLOTS_XFLR5;
Plot_Options.PRINT_PLOTS_AERO = PRINT_PLOTS_AERO;
Plot_Options.mark_Type = mark_Type;
Plot_Options.mark_legend = mark_legend;
Plot_Options.VECTOR = VECTOR;
Plot_Options.VECTOR_XFLR5 = VECTOR_XFLR5;
PL = 1;

if PRINT_PLOTS_XFLR5 == 1
    [Fig] = plot_XFLR5_Aero_Aug2018(Aero,DATA_Ae,DATA_PL,Plot_Options,Fig,Design_criteria);
    [Fig] = plot_XFLR5_Polar_Aug2018(Aero,DATA_Ae,DATA_PL,Plot_Options,Fig,Design_criteria)
    pause
elseif PRINT_PLOTS_AERO == 1
    [Fig] = plot_Aero_Aug2018(Aero,DATA_Ae,DATA_PL,V_Performance,Plot_Options,Fig);
    pause
else
%     Fig = 0;
end

% % Actualze Geometry 
% % Correction of location of mean aerodynamic center (MAC) using XFLR5 
% MAC.xbar_w1 = 0.057;
% MAC.xbar_w2 = 0.056;
% MAC.xbar_v = 0.083;
% 
% Geo_tier = Geometric_Data_tier_UPDATE(Geo_tier,MAC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates Geometry Geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XFLR5_file = strcat('fuse_pepi_desplazado.txt');
XFLR5_file = strcat('Fuselage_DATA_scaled.txt');

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
[Body_Geo,meshData] = Calc_Body_Geometry_Dec2018(XFLR5_file,Geo_tier,nPoints,nSections,conversion_UNITS,lecture_Geo);

% Fusion Aerodynamic Properties
[Aero_TH] = Polar_Literatura(h,V,Geo_tier,conv_UNITS,Body_Geo);


% Generates Data for Ploting 3D Mesh
% [Fig] = plot_GEOMETRY_2018(Geo_tier,Plot_Options,Body_Geo,meshData,Propulsion);

% DynVar = Get_ForcesMoments(x,delta,Data_Trim,Data_Der,h,V)

% DATA_int,Data_Der,Ctrl_Sig,Data_ATM,Wind,Data_ref)

