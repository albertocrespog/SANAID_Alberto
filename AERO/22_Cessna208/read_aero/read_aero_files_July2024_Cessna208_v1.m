function [casos prefix mark_legend X_OC] = read_aero_files_July2024_Cessna208_v1(Performance)

V = Performance.V;
h = Performance.h;
 
% Wing 1 Cessna 208
casE = 1;
% wing
casos{casE} = 'Wing_Cessna208_T1-90_0 m_s-LLT-x0_56m.txt'; 
X_OC(casE) = 0.56;
st{casE} = strcat('LLT - W_1');

casE = casE + 1;
casos{casE} = 'Wing_Cessna208_T1-90_1 m_s-VLM1-x0_56m_FL100_wingtip_pressure.txt'; 
X_OC(casE) = 0.56;
st{casE} = strcat('VLM1 - W_1');

casE = casE + 1;
casos{casE} = 'Wing_Cessna208_T1-90_1 m_s-TriUniform-ThickSurf-x0_58m_FL100.txt'; 
X_OC(casE) = 0.58;
st{casE} = strcat('TriUniform - W_1');


% HTP
casE = casE + 1;
% HTP & VTP
casos{casE} = 'HTPVTP_T1-90_1 m_s-LLT-x0_39m_FL100.txt'; % Case 1
X_OC(casE) = 0.39;
st{casE} = strcat('TriUniform - HTP_1');

casE = casE + 1;
% HTP & VTP
casos{casE} = 'HTPVTP_T1-90_1 m_s-VLM1-x0_3885m_FL100.txt'; % Case 1
X_OC(casE) = 0.3885;
st{casE} = strcat('VLM1 - HTP_1');

casE = casE + 1;
% HTP & VTP
casos{casE} = 'HTPVTP_T1-90_1 m_s-TriUniform-ThinSurf-x0_3825m_FL100.txt'; % Case 1
X_OC(casE) = 0.3825;
st{casE} = strcat('TriUniform - HTP_1');

% VTP lateral
casE = casE + 1;
casos{casE} = 'VTP_sidelslip_T5-a0_0°-90_1m_s-VLM2-x0_4m.txt'; % Case 2
X_OC(casE) = 0.00;
st{casE} = strcat('VLM2 - VTP');

prefix = strcat('Cessna208_'); 

mark_legend = st;


