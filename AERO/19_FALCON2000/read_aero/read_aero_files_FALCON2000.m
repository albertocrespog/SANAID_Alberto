function [casos prefix mark_legend X_OC] = read_aero_files_FALCON2000(Performance)

V = Performance.V;
h = Performance.h;
 
% Wing Complete alone
casE = 1;
% wing
casos{casE} = 'T1-20_0 m_s-LLT-x0_52m_3000m.txt'; % Case 1
X_OC(casE) = 0.085;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

% VTP Longitudinal
casE = casE + 1;
casos{casE} = 'T5-a0_0-20_0m_s-VLM2-x0_5m.txt'; % Case 2
X_OC(casE) = 0.06475;
st{casE} = strcat('VLM - VTP (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

prefix = strcat('SOLARTII_'); 

mark_legend = st;