function [casos prefix mark_legend X_OC] = read_aero_files_Dec2024_H2_v1(Performance)

V = Performance.V;
h = Performance.h;
 
% Wing
casE = 1;
% wing
casos{casE} = 'wing_T1-30_0 m_s-TriUniform-ThickSurf-x0_04m.txt'; 
X_OC(casE) = 0.04;
st{casE} = strcat('TriUniform - W_1');

% wing - lateral
casE = casE + 1;
casos{casE} = 'wing_T5-a0_0°-30_0m_s-TriUniform-ThickSurf-x0_04m.txt'; 
X_OC(casE) = 0.04;
st{casE} = strcat('TriUniform - W_1');

% Wing 2
casE = casE + 1;
% wing
casos{casE} = 'can_T1-30_0 m_s-TriUniform-ThickSurf-x0_04m.txt'; 
X_OC(casE) = 0.04;
st{casE} = strcat('TriUniform - W_1');

% wing 2 - lateral
casE = casE + 1;
casos{casE} = 'can_T5-a0_0°-30_0m_s-TriUniform-ThickSurf-x0_04m.txt'; 
X_OC(casE) = 0.04;
st{casE} = strcat('TriUniform - W_1');

% VTP
casE = casE + 1;
% wing
casos{casE} = 'w1VTP_T1-30_0 m_s-TriUniform-ThickSurf-x0_04m.txt'; % Case 1
X_OC(casE) = 0.0475;
st{casE} = strcat('TriUniform - HTP_1');

% VTP lateral
casE = casE + 1;
casos{casE} = 'w1VTP_T5-a0_0°-30_0m_s-TriUniform-ThickSurf-x0_04m.txt'; % Case 2
X_OC(casE) = 0.0475;
st{casE} = strcat('TriUniform - VTP');

prefix = strcat('H2_'); 

mark_legend = st;


