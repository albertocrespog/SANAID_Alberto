function [casos prefix mark_legend X_OC] = read_aero_files_MILVUS(Performance)

V = Performance.V;
h = Performance.h;

% Estudios
% Front Wing analysis for different locations
casE = 1;
casos{casE} = 'MILVUS_T1-25_0_m_s-LLT-x0_02875m.txt'; 
X_OC(casE) = 0.02875;
st{casE} = strcat('LLT - can (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
casE = casE + 1;
% Rear Wing analysis for different locations (same as front)
casos{casE} = 'MILVUS_T1-25_0_m_s-LLT-x0_02875m.txt'; 
X_OC(casE) = 0.02875;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
casE = casE + 1;
% Vtail
casos{casE} = 'MILVUS_T1-25_0_m_s-VLM2_Vee.txt'; 
X_OC(casE) = 0.000;
st{casE} = strcat('VLM2 - Vee (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
casE = casE + 1;

% Sideslip
casos{casE} = 'MILVUS_T5-a0_0-25_0m_s-VLM2-x0_02875m.txt'; 
X_OC(casE) = 0.000;
st{casE} = strcat('VLM2 - can (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
casE = casE + 1;

casos{casE} = 'MILVUS_T5-a0_0-25_0m_s-VLM2-x0_02875m.txt'; 
X_OC(casE) = 0.000;
st{casE} = strcat('VLM2 - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
casE = casE + 1;

casos{casE} = 'MILVUS_T5-a0_0-25_0m_s-VLM2_Vee.txt'; 
X_OC(casE) = 0.000;
st{casE} = strcat('VLM2 - Vee (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
casE = casE + 1;

prefix = strcat('MILVUS_');
mark_legend = st;