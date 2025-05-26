function [casos prefix mark_legend X_OC] = read_aero_files_A400(Performance)

V = Performance.V;
h = Performance.h;
 
% Wing Complete alone
casE = 1;
% wing
casos{casE} = 'T6_1-WING-TriUniform-ThickSurf-Vel(200_65)_x_4.8m.txt'; % Case 1
X_OC(casE) = 4.8;
st{casE} = strcat('TriUniform - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

% HTP & VTP Longitudinal
casE = casE + 1;
casos{casE} = 'T6_1-HTP&VTP-TriUniform-ThickSurf-Vel(200_65,)-alpha_x_2.885m.txt'; % Case 2
X_OC(casE) = 2.885;
st{casE} = strcat('TriUniform - HTP (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

% HTP & VTP Lateral-directional
casE = casE + 1;
casos{casE} = 'T6_1-HTP&VTP-VLM2-Vel(200_65)-beta-CGx(2_88).txt';  % Case 3
X_OC(casE) = 2.88;
st{casE} = strcat('VLM2 - VTP (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

prefix = strcat('A400_'); 

mark_legend = st;