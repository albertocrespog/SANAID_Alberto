function [casos prefix mark_legend X_OC] = read_aero_files_Nov2024_HA_v1(V_Performance)

V = V_Performance.V;
h = V_Performance.h;

% Estudios
% 1 - Front Wing analysis for different locations
casE = 1;
casos{casE} = '26_HA_wing_T1-35_0 m_s-TriUniform-ThickSurf-x0_0675m.txt';
X_OC(casE) = 0.0675;
st{casE} = strcat('TRIUNIFORM - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');

% 2 - Front Wing analysis for different locations
% Wing - Sideslip
casE = casE + 1;
casos{casE} = '26_HA_wing_T5-a0_0°-35_0m_s-VLM2-x0_0675m.txt'; 
X_OC(casE) = 0.0675;
st{casE} = strcat('VLM2 - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');

% 3 - Vtail UP - Alone
casE = casE + 1;
casos{casE} = '26_HA_vee_T1-35_0 m_s-TriUniform-ThinSurf-x0_09m-z0_15m.txt';
X_OC(casE) = 0.090;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');

% 4 - Vtail UP sideslip
casE = casE + 1;
casos{casE} = '26_HA_vee_T5-a0_0°-35_0m_s-TriUniform-ThinSurf-x0_09m-z0_15m.txt'; % Flap
X_OC(casE) = 0.090;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');


prefix = strcat('HA_');
mark_legend = st;