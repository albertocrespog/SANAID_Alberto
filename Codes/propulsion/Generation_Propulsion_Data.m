function [Prop_data] = Generation_Propulsion_Data(AC_CONFIGURATION,OUTPUT_read_XLSX,filenameS)

SF = OUTPUT_read_XLSX.AC_Data_flags.SF;

%% Propulsive Model
%% Compares 3 prop models
% - Model 1 - APC data
% - Model 2 - Wind tunnel data for different props
% 1 - APC 20x8
% 2 - APC 22x10
% 3 - APC 22x12
% 4 - APC 22x12W % 22x12WE Thin Electric Wide Blade - https://www.apcprop.com/product/20x12we/
% 5 - APC 21x13W
% 6 - APC 21x14
% - Model 3 - Wind tunnel data for different angle of attack for APC 22x12W
model_prop = OUTPUT_read_XLSX.Propulsive_flags.model_prop; %
prop_selec_APC = OUTPUT_read_XLSX.Propulsive_flags.prop_selec_APC; %
prop_selec_WT1 = OUTPUT_read_XLSX.Propulsive_flags.prop_selec_WT1; %
prop_selec_WT2 = OUTPUT_read_XLSX.Propulsive_flags.prop_selec_WT2; %
Prop_selection.model_prop = model_prop;
% Selects the prop for
% - Model 1 - APC data
% - Model 2 - Wind tunnel data for different props
% Stores info
Prop_selection.prop_selec_APC = prop_selec_APC;
Prop_selection.prop_selec_WT1 = prop_selec_WT1;
Prop_selection.prop_selec_WT2 = prop_selec_WT2;


%% Compares 3 prop models
% - Model 1 - APC data 
% - Model 2 - Wind tunnel data for different props
% - Model 3 - Wind tunnel data for different angle of attack for APC 22x12W
model_prop = Prop_selection.model_prop;
% Selects the prop for 
% - Model 1 - APC data 
% - Model 2 - Wind tunnel data for different props
prop_selec_APC = Prop_selection.prop_selec_APC;
prop_selec_WT1 = Prop_selection.prop_selec_WT1;

% Propeller Diameter
D_prop = OUTPUT_read_XLSX.Propulsive_flags.D_prop; %

%% Engine Configuration
n_eng = AC_CONFIGURATION.n_eng;
Prop_type = AC_CONFIGURATION.Prop_type;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Engine properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Propulsive model
% AXI 5360/24HD V2 GOLD LINE
l_eng = SF*0.104; %length
d_eng = SF*0.063; %diameter

Prop_data.l_eng = l_eng;
Prop_data.d_eng = d_eng;
Prop_data.n_eng = n_eng;

% Approximate dimenstions of the nacelles
l_nc = l_eng*2;
d_nc= d_eng*2; 

% e = Engine number concerning to excel master table. From 1 to 6.
e=6;
load('Data.mat')
load('Data_Prop.mat')

%---------------------------ENGINE DATA----------------------------------%
W_eng = Data.engine(e).Weight;     % Engine Weight.
D_propin=D_prop*100/2.54; % Propeller Diameter [in].
A_prop=pi*(D_prop/2)^2;          % Swept area by propeller.

% Propeller data
% Datos genéricos Hélice para 22x12W - They are used to scale the Prop
b_p = (22*2.54/100);
c_p = 3/100;
b_p_c_p = b_p/c_p;
c_prop = D_prop/b_p_c_p;
S_prop = D_prop*c_prop;
AR_prop = (D_prop^2)/S_prop;
RPM_max = 150000/D_propin; % Max RPM by engine builder - https://www.apcprop.com/technical-information/rpm-limits/
% Thin Electric (E) Propellers
% Maximum RPM=150000/prop diameter (inches)
Prop_data.W_eng = W_eng;
Prop_data.D_prop = D_prop;
Prop_data.D_propin = D_propin;
Prop_data.A_prop = A_prop;
Prop_data.RPM_max = RPM_max;
Prop_data.c_prop = c_prop;
Prop_data.S_prop = S_prop;
Prop_data.AR_prop = AR_prop;

Prop_data.l_eng = l_eng;
Prop_data.d_eng = d_eng;
Prop_data.l_nc = l_nc;
Prop_data.d_nc = d_nc;

% Efficiencies
eta_gear=0.96;                   % Gear box efficiency.
eta_m=0.88;                      % Engine efficiency (output/input).
eta_esc=0.98;                    % Speed controller efficiency.
eta_dist=0.96;                   % Shaft efficiency.
Prop_data.eta_gear = eta_gear;
Prop_data.eta_m = eta_m;
Prop_data.eta_esc = eta_esc;
Prop_data.eta_dist = eta_dist;

%% Propulsive Model
%% Compares 3 prop models
% - Model 1 - APC data 
% - Model 2 - Wind tunnel data for different props
% 1 - APC 20x8
% 2 - APC 22x10 
% 3 - APC 22x12 
% 4 - APC 22x12W % 22x12WE Thin Electric Wide Blade - https://www.apcprop.com/product/20x12we/
% 5 - APC 21x13W 
% 6 - APC 21x14  
% - Model 2 - Wind tunnel data for different angle of attack for APC 22x12W

%% Propulsive Model - Model 1 - APC data 
% Reads APC data file
[casos_prop_APC prefix_prop mark_legend_prop]= read_prop_files_May2020;
[Prop_data_APC] = Read_data_prop_APC(casos_prop_APC);
PROPDATA_APC = process_files_prop_APC(casos_prop_APC,Prop_data_APC);

%% Propulsive Model - Model 2 - Wind tunnel data for different props
% 1 - APC 20x8
% 2 - APC 22x10 
% 3 - APC 22x12 
% 4 - APC 22x12W % 22x12WE Thin Electric Wide Blade - https://www.apcprop.com/product/20x12we/
% 5 - APC 21x13W 
% 6 - APC 21x14  
N_props_WT1= 6;
[casos_prop_WT1 prefix mark_legend] = read_prop_files_WT1;
for i=1:N_props_WT1
    Prop_type = i; 
    Prop_data_WT1{i} = Prop_data_Wind_tunnel_1(Prop_type);
    PROPDATA_WT1{i} = process_files_prop_WT1(casos_prop_WT1,Prop_data_WT1{i},N_props_WT1);
end

%% Propulsive Model - Model 2 - Wind tunnel data for different props
% APC 22x12W % 22x12WE Thin Electric Wide Blade - https://www.apcprop.com/product/20x12we/
[Prop_data_WT2] = Read_data_prop_WT2;
PROPDATA_WT2 = process_files_prop_WT2(Prop_data_WT2);

prop_selec_APC = Prop_selection.prop_selec_APC;
prop_selec_WT1 = Prop_selection.prop_selec_WT1;

switch model_prop
    case 1 % - Model 1 - APC data
        k = prop_selec_APC;
        Prop_data.CT_Polyfit = PROPDATA_APC.CT_Polyfit_APC{k};
        Prop_data.CP_Polyfit = PROPDATA_APC.CP_Polyfit_APC{k};
        Prop_data.CQ_Polyfit = PROPDATA_APC.CQ_Polyfit_APC{k};
        Prop_data.etamp_Polyfit = PROPDATA_APC.etap_Polyfit_APC{k};
        Prop_data.N_order_CT = PROPDATA_APC.N_order_CT{k};
        Prop_data.N_order_CP = PROPDATA_APC.N_order_CP{k};
        Prop_data.N_order_CQ = PROPDATA_APC.N_order_CQ{k};
        Prop_data.N_order_etamp = PROPDATA_APC.N_order_etap{k};
        
    case 2  % - Model 2 - Wind tunnel data for different props
        k = prop_selec_WT1;
        Prop_data.CT_Polyfit = PROPDATA_WT1{k}.CT_Polyfit_WT1;
        Prop_data.CP_Polyfit = PROPDATA_WT1{k}.CP_Polyfit_WT1;
        Prop_data.CQ_Polyfit = PROPDATA_WT1{k}.CQ_Polyfit_WT1;
        Prop_data.etamp_Polyfit = PROPDATA_WT1{k}.etamp_Polyfit_WT1;
        Prop_data.N_order_CT = PROPDATA_WT1{k}.N_order_CT;
        Prop_data.N_order_CP = PROPDATA_WT1{k}.N_order_CP;
        Prop_data.N_order_CQ = PROPDATA_WT1{k}.N_order_CQ;
        Prop_data.N_order_etamp = PROPDATA_WT1{k}.N_order_etamp;
        
    case 3 % - Model 3 - Wind tunnel data for different angle of attack for APC 22x12W
        % Future versions
end

Saving_data_Prop(Prop_data,OUTPUT_read_XLSX,filenameS)
% prefixa = strcat(OUTPUT_read_XLSX.PLOT_flags.fname);
% st0 = strcat('data\');
% st1 = strcat(st0,prefixa);
% st2 = strcat('\Prop_data.mat');
% name   = strcat(st1,st2);
% save(name, 'Prop_data')
% save('data/Prop_data.mat', 'Prop_data')
