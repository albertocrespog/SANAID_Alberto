function [rho_fus_fairing,rho_w,rho_can,rho_HTP,rho_VTP,rho_Vee] = EMERGENTIadensities(AC_CONFIGURATION,OUTPUT_read_XLSX,Body_Geo,Geo_tier); % calculated from ProVANT-EMERGENTIA 100% I

%% identifies the aerodynamic surfaces being used
W1 = AC_CONFIGURATION.W1;
HTP = AC_CONFIGURATION.HTP;
VTP = AC_CONFIGURATION.VTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
case_AC = AC_CONFIGURATION.case_AC;

% Weight of structure EMRGENTIA
m_w1 = OUTPUT_read_XLSX.Weights_flags.m_w1;
m_HTP = OUTPUT_read_XLSX.Weights_flags.m_HTP;
m_VTP = OUTPUT_read_XLSX.Weights_flags.m_VTP;
m_Can = OUTPUT_read_XLSX.Weights_flags.m_Can;
m_Vee = OUTPUT_read_XLSX.Weights_flags.m_Vee;
m_Vee2 = OUTPUT_read_XLSX.Weights_flags.m_Vee2;
m_fus_fairing = OUTPUT_read_XLSX.Weights_flags.m_fus_fairing;

% Fudge factor
f_f = OUTPUT_read_XLSX.Weights_flags.f_f;
% Scaling Factor
SF = OUTPUT_read_XLSX.AC_Data_flags.SF;

if W1 == 1
    S_w1_s = Geo_tier.S_w1_s;
    rho_w =m_w1/S_w1_s; 
    m_w1 = f_f*rho_w*S_w1_s;                             % Wing weight.
else
    m_w1 = 0;
    rho_w = 0
end
if HTP == 1
    S_HTP_s = Geo_tier.S_HTP_s;
    rho_HTP =m_HTP/S_HTP_s; 
    m_HTP = f_f*rho_HTP*S_HTP_s;                    %Peso HTP
else
    m_HTP = 0;
    rho_HTP = 0;
end
if VTP == 1
    S_VTP_s = Geo_tier.S_VTP_s;
    rho_VTP =m_VTP/S_VTP_s; 
    m_VTP = f_f*rho_VTP*S_VTP_s;                    %Peso HTP
else
    m_VTP = 0;
    rho_VTP = 0;
end
if Can == 1
    S_can_s = Geo_tier.S_can_s;
    rho_can = rho_w; 
    m_Can = f_f*rho_can*S_can_s;                  %Peso canard
else
    m_Can = 0;
    rho_can = 0;
end
if Vee == 1
    S_vee_s = Geo_tier.S_vee_s;
    rho_Vee = m_Vee/S_vee_s; 
    m_Vee = f_f*rho_Vee*S_vee_s;                 %Peso cola en v
else
    m_Vee = 0;
    rho_Vee = 0;
end
if Vee2 == 1
    S_vee_s = Geo_tier.S_vee_s;
    rho_Vee = m_Vee/S_vee_s; 
    m_Vee = f_f*rho_Vee*S_vee_s;                 %Peso cola en v
else
    m_Vee = 0;
    rho_Vee = 0;
end

% Data fron Céfiro III for carbon fiber composite fairing 
Surf_TOT = Body_Geo.Surf_TOT; % m^2

% Recalculo de las nuevas densidades
rho_fus_fairing = m_fus_fairing/Surf_TOT;