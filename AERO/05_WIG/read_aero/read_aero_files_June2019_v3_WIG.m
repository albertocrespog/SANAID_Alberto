function [casos prefix mark_legend X_OC] = read_aero_files_Feb2020_v3_WIG(Performance)

V = Performance.V;
h = Performance.h;

% Estudios
% Front Wing analysis for different locations
casE = 1;
casos{casE} = 'W1_T1-125_0x2_3m-TG_GE5.txt'; 
X_OC(casE) = 2.3;
st{casE} = strcat('VLM - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
casE = casE + 1;
casos{casE} = 'HTP&VTP_T1-125_0x0_401m-TG_GE9.txt'; 
X_OC(casE) = 0.401;
st{casE} = strcat('VLM - HTP & VTP (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
casE = casE + 1;
casos{casE} = 'HTP&VTP_betaT5-125_0_x0_401m-TG_GE9.txt'; 
X_OC(casE) = 0.401;
st{casE} = strcat('VLM - HTP & VTP @ sideslip (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
casE = casE + 1;

prefix = strcat('WIG_');

mark_legend = st;