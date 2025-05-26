function [casos prefix mark_legend X_OC] = read_aero_files_ALO_v1(Performance)

V = Performance.V;
h = Performance.h;

% Estudios
% Front Wing analysis for different locations
casE = 1;
casos{casE} = 'T1-32_0 m_s-VLM2-3_0kg-x0_221mALA.txt'; 
X_OC(casE) = 0.000;
st{casE} = strcat('VLM - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
casE = casE + 1;
casos{casE} = 'T1-32_0 m_s-VLM2-0_8kg-x0_096m-z0_1mELEVATOR.txt'; 
X_OC(casE) = 0.005;
st{casE} = strcat('VLM - W_2 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
casE = casE + 1;

prefix = strcat('ALO_');


mark_legend = st;