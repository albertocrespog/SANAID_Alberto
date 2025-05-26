function Design_criteria = Generate_aerodynamic_data(OUTPUT_read_XLSX,AC_CONFIGURATION,conv_UNITS) 

R2D = conv_UNITS.R2D;
D2R = conv_UNITS.D2R;

%% identifies the aerodynamic surfaces being used
W1 = AC_CONFIGURATION.W1;
HTP = AC_CONFIGURATION.HTP;
VTP = AC_CONFIGURATION.VTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
% Nac = AC_CONFIGURATION.Nac;

%% Generates the file for Aerodynamic Data
index_w1 = OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1;
index_can = OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can;
index_vee = OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee;
index_vee2 = OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee2;
index_HTP = OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_HTP;
index_VTP = OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_VTP;
%
index_w1_cy = OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1_cy; %
index_HTP_cy = OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_HTP_cy; %
index_vee_cy = OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_cy; %
index_vee2_cy = OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee2_cy; %
index_can_cy = OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can_cy; %
index_VTP_cy = OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_VTP_cy; %
%
index_w1_VLM = OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_w1_VLM; %
index_HTP_VLM = OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_HTP_VLM; %
index_vee_VLM = OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee_VLM; %
index_vee2_VLM = OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_vee2_VLM; %
index_can_VLM = OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_can_VLM; %
index_VTP_VLM = OUTPUT_read_XLSX.Aerodynamic_Data_flags.index_VTP_VLM; %
%
i_w1 = OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_w1; %
i_HTP = OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_HTP; %
i_vee = OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_vee; %
i_vee2 = OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_vee2; %
i_can = OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_can; %
i_VTP = OUTPUT_read_XLSX.Aerodynamic_Data_flags.i_VTP; %
%
       
if W1 == 1
    %% Stores data
    alpha_selected_w1 = i_w1; % (degs)
    Design_criteria.index_w1 = index_w1;
    Design_criteria.index_w1_VLM = index_w1_VLM;
    Design_criteria.index_w1_cy = index_w1_cy;
    Design_criteria.alpha_selected_w1 = alpha_selected_w1;
    % Design choices for the incidence of the different surfaces
    Design_criteria.i_w1 = i_w1*D2R; % incidence of Front Wing
end

if Can  == 1
    %% Stores data
    alpha_selected_can = i_can; % (degs)
    Design_criteria.index_can = index_can;
    Design_criteria.index_can_VLM = index_can_VLM;
    Design_criteria.index_can_cy = index_can_cy;
    Design_criteria.alpha_selected_can = alpha_selected_can;
    % Design choices for the incidence of the different surfaces
    Design_criteria.i_can = i_can*D2R; % incidence of Front Wing
end

if HTP == 1
    %% Stores data
    alpha_selected_HTP = i_HTP; % (degs)
    Design_criteria.index_HTP = index_HTP;
    Design_criteria.index_HTP_VLM = index_HTP_VLM;
    Design_criteria.index_HTP_cy = index_HTP_cy;
    Design_criteria.alpha_selected_HTP = alpha_selected_HTP;
    % Design choices for the incidence of the different surfaces
    Design_criteria.i_HTP = i_HTP*D2R; % incidence of Rear Wing
end

if Vee == 1
    %% Stores data
    alpha_selected_vee = i_vee; % (degs)
    Design_criteria.index_vee = index_vee;
    Design_criteria.index_vee_VLM = index_vee_VLM;
    Design_criteria.index_vee_cy = index_vee_cy;
    Design_criteria.alpha_selected_vee = alpha_selected_vee;
    % Design choices for the incidence of the different surfaces
    Design_criteria.i_vee = i_vee*D2R; % incidence of Rear Wing
end

if Vee2 == 1
    %% Stores data
    alpha_selected_vee2 = i_vee2; % (degs)
    Design_criteria.index_vee2 = index_vee2;
    Design_criteria.index_vee2_cy = index_vee2_cy;
    Design_criteria.index_vee2_VLM = index_vee2_VLM;
    Design_criteria.alpha_selected_vee2 = alpha_selected_vee2;
    % Design choices for the incidence of the different surfaces
    Design_criteria.i_vee2 = i_vee2*D2R; % incidence of Rear Wing
end

if VTP == 1
   alpha_selected_VTP = i_VTP; % (degs)
    Design_criteria.index_VTP = index_VTP;
    Design_criteria.index_VTP_VLM = index_VTP_VLM;
    Design_criteria.index_VTP_cy = index_VTP_cy;
    Design_criteria.alpha_selected_VTP = alpha_selected_VTP;
    % Design choices for the incidence of the different surfaces
    Design_criteria.i_VTP = i_VTP*D2R; % incidence of Rear Wing
end

%% Flight Safe Margin to calculate aerodynamic properties
Flight_SF = OUTPUT_read_XLSX.Performance_pre_flags.Flight_SF;
Design_criteria.Flight_SF = Flight_SF;
%% Revisar si pasar al Excel

% %% Propulsion Generation
% alpha_f = 0*D2R;
% beta_f = 0*D2R;