function [casos prefix mark_legend X_OC] = read_aero_files_July2024_Future4_v1_FR1(Performance)

V = Performance.V;
h = Performance.h;
 
% Wing simple 
casE = 1;
% wing
casos{casE} = 'w1_T1-200_0 m_s-TriUniform-ThickSurf-x0_095m.txt'; % Case 1
X_OC(casE) = 0.095;
st{casE} = strcat('TriUniform - W_1');

% Wing lateral
casE = casE + 1;
casos{casE} = 'w1_T5-a0_0°-200_0m_s-VLM2-x0_095m.txt'; % Case 3
X_OC(casE) = 0.095;
st{casE} = strcat('VLM2 - W_1');


% VTail
casE = casE + 1;
casos{casE} = 'Vee_T1-200_0 m_s-TriUniform-ThickSurf-x0_085m-z0_05m.txt'; % Case 2
X_OC(casE) = 0.085;
st{casE} = strcat('TriUniform - Vee');

% VTail lateral
casE = casE + 1;
casos{casE} = 'Vee_T5-a0_0°-200_0m_s-TriUniform-ThickSurf-x0_085m-z0_05m.txt'; % Case 4
X_OC(casE) = 0.085;
st{casE} = strcat('VLM2 - Vee');

prefix = strcat('FR1_'); 

mark_legend = st;


