function [casos prefix mark_legend X_OC] = read_aero_files_July2024_KingAir350_v1(Performance)

V = Performance.V;
h = Performance.h;
 
% Wing 1 Cessna 208
casE = 1;
% wing
casos{casE} = 'w1_T1-100_3 m_s-TriUniform-ThickSurf-x0_8m.txt'; 
X_OC(casE) = 0.8;
st{casE} = strcat('TriUniform - W_1');

casE = casE + 1;
casos{casE} = 'w1_T1-84_9 m_s-TriUniform-ThickSurf-x0_8m.txt'; 
X_OC(casE) = 0.8;
st{casE} = strcat('TriUniform - W_1');

casE = casE + 1;
casos{casE} = 'w1_T1-84_9 m_s-TriUniform-ThickSurf-x0_8m_SL.txt'; 
X_OC(casE) = 0.8;
st{casE} = strcat('TriUniform - W_1');

casE = casE + 1;
casos{casE} = 'w1_T5-a0_0-100_3m_s-TriUniform-ThickSurf-x0_8m.txt'; 
X_OC(casE) = 0.8;
st{casE} = strcat('TriUniform - W_1');

% HTP
casE = casE + 1;
% HTP
casos{casE} = 'HTP_T1-100_3 m_s-TriUniform-ThickSurf-x0_8m.txt'; % Case 1
X_OC(casE) = 0.8;
st{casE} = strcat('TriUniform - HTP_1');

casE = casE + 1;
% HTP
casos{casE} = 'HTP_T1-84_9 m_s-TriUniform-ThickSurf-x0_8m_SL.txt'; % Case 1
X_OC(casE) = 0.3885;
st{casE} = strcat('VLM1 - HTP_1');

casE = casE + 1;
% HTP
casos{casE} = 'HTP_T5-a0_0-84_9m_s-TriUniform-ThickSurf-x0_8m_SL.txt'; % Case 1
X_OC(casE) = 0.3825;
st{casE} = strcat('TriUniform - HTP_1');

% VTP lateral
casE = casE + 1;
casos{casE} = 'HTP_VTP_T1-100_3 m_s-TriUniform-ThickSurf-x0_8m.txt'; % Case 2
X_OC(casE) = 0.08;
st{casE} = strcat('TriUniform - VTP');

casE = casE + 1;
casos{casE} = 'VTP_T5-a0_0-100_3m_s-TriUniform-ThickSurf-x0_8m.txt'; % Case 2
X_OC(casE) = 0.08;
st{casE} = strcat('TriUniform - VTP');



prefix = strcat('KingAir_'); 

mark_legend = st;
