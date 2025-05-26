close all
clear all
% Defines Generic Path so that files can be organiez in folders
% change for the different versions
addpath(genpath('../EMERGENTIA'))
% Units conversion
conv_UNITS = conversion_UNITS;

% Flight Conditons
h = 200;
V = 25;
[Temp,rho,p,a]=atmos_inter(h);

Data_ATM.Temp = Temp;
Data_ATM.rho = rho;
Data_ATM.p = p;
Data_ATM.a = a;

% Obtains info from Cefiro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NEED TO BE REPLACED BY EMERGENTIA
%
Data_Trim = Trim_Cefiro_1214(h,V);
Data_Der = get_derivatives_Cefiro_1214(Data_Trim,h,V);

% Routine that generates States
[x,delta]=Generate_States(Data_Trim,Data_Der,h,V);
D = 0.6043;
n = 5000/60;
alpha = Data_Trim.alpha_1;
% Number of prop
k = 4;
% Plots figures for Engine model
Plot_Options.figures_CT_model = 0;
fig = 1;
% Gets the engine properties
[Propulsion] = get_EngineProperties(Data_Trim,V,alpha,n,D,k,fig,Plot_Options); 
% Initial Geometry
fig = 0;
Geo_tier = Geometric_Data_tier;

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

T = Propulsion.Ti;
D_prop = 0.6043;
Sheli=(pi*(D_prop/2)^2);              %superficie de la hélice
vio=sqrt(T/(2*rho*Sheli));          %velocidad inicial inducida por hélice
vi = -0.5*V + sqrt(0.25*V^2 + vio^2);    %velocidad real inducida por la hélice
V_nc = V + vi     ;                      %velocidad que afecta a la cola
q_prop=0.5*rho*V_nc^2;  
% Drag coefficient nacelle
CD0_nc = Aero_TH.CD0_nc;
% Drag nacelle
D_nc =q_prop*Geo_tier.S_w1*CD0_nc;

% Wind components
Wx = 0;
Wy = 0;
Wz = 0;

Wind.Wx = Wx;
Wind.Wy = Wy;
Wind.Wz = Wz;
DynVar = Get_ForcesMoments(x,delta,Data_Trim,Data_Der,h,V,Wind,Data_ATM,Propulsion)
% DATA_int,Data_Der,Ctrl_Sig,Data_ATM,Wind,Data_ref),

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stores information for the plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data for the legend
% V_CR = V_Performance.V_CR; % Cruise speed
% h_i = V_Performance.h_i; % Altitude
% Generates the Legend
[mark_Type,mark_legend] = Generates_plot_Info(V,h); % Generates the style of the lines and the Legend

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
Plot_Options.mark_Type = mark_Type;
Plot_Options.mark_legend = mark_legend;
Plot_Options.VECTOR = VECTOR;
Plot_Options.VECTOR_XFLR5 = VECTOR_XFLR5;
PL = 1;

% Generates Data for Ploting 3D Mesh
[Fig] = plot_GEOMETRY_2018(Geo_tier,Plot_Options,Body_Geo,meshData,Propulsion);


