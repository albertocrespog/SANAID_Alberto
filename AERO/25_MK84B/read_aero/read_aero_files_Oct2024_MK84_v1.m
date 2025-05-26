function [casos prefix mark_legend X_OC] = read_aero_files_Oct2024_MK84_v1(V_Performance)

V = V_Performance.V;
h = V_Performance.h;

% Estudios
% 1 - Front Wing analysis for different locations
casE = 1;
casos{casE} = 'MK84_T1-200_0 m_s-TriUniform-ThickSurf-x0_181m.txt';
X_OC(casE) = 0.181;
st{casE} = strcat('VLM2- W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');

% 2 - Front Wing analysis for different locations
% Wing - Sideslip
casE = casE + 1;
casos{casE} = 'MK84_T5-a0_0°-200_0m_s-VLM2-x0_181m.txt'; 
X_OC(casE) = 0.181;
st{casE} = strcat('LLT - W_2 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');

% 3 - Vtail UP - Alone
casE = casE + 1;
casos{casE} = 'MK84_Vtail_DOWN_T1-200_0 m_s-TriUniform-ThickSurf-x0_0625m-z-0_212m.txt';
X_OC(casE) = 0.212;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');

% 4 - Vtail UP sideslip
casE = casE + 1;
casos{casE} = 'MK84_Vtail_DOWN_T5-a0_0°-200_0m_s-TriUniform-ThickSurf-x0_0625m-z-0_212m.txt'; % Flap
X_OC(casE) = 0.212;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');

% 3 - Vtail inverted - Alone
casE = casE + 1;
casos{casE} = 'MK84_Vtail_UP_T1-200_0 m_s-TriUniform-ThickSurf-x0_212m.txt';
X_OC(casE) = 0.212;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');

% 4 - Vtail inverted sideslip
casE = casE + 1;
casos{casE} = 'MK84_Vtail_UP_T5-a0_0°-200_0m_s-TriUniform-ThickSurf-x0_212m.txt'; % Flap
X_OC(casE) = 0.212;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');


prefix = strcat('MK84_');
mark_legend = st;