function [casos prefix mark_legend X_OC] = read_aero_files_July2024_Future7_v1(Performance)

V = Performance.V;
h = Performance.h;
 
% Wing 1 CIRIO
casE = 1;
% wing
casos{casE} = 'T1-24_0 m_s-TriUniform-ThinSurf-x0_07m.txt'; 
X_OC(casE) = 0.07;
st{casE} = strcat('TriUniform - W_1');

% HTP
casE = casE + 1;
% wing
casos{casE} = 'T1-24_0 m_s-TriUniform-ThinSurf-x0_73_HTP.txt'; % Case 1
X_OC(casE) = 0.073;
st{casE} = strcat('TriUniform - HTP_1');

% VTP lateral
casE = casE + 1;
casos{casE} = 'T5-a0_0-24_0m_s-TriUniform-ThinSurf-x0_0m_VTP.txt'; % Case 2
X_OC(casE) = 0.00;
st{casE} = strcat('TriUniform - VTP');

prefix = strcat('SOLARTII_'); 

mark_legend = st;


