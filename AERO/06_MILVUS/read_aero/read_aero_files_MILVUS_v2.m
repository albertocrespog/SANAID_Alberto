function [casos prefix mark_legend X_OC] = read_aero_files_MILVUS_v2(Performance)

V = Performance.V;
h = Performance.h;

% Estudios
% Front Wing analysis for different locations
casE = 1;
casos{casE} = 'wing_1-25_0 m_s-TriUniform-ThickSurf-x0_03m.txt'; 
X_OC(casE) = 0.03;
st{casE} = strcat('TRI - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
casE = casE + 1;

% Front Wing analysis SideSlip (same as front)
casos{casE} = 'wing_T5-a0_0°-25_0m_s-TriUniform-ThickSurf-x0_03m.txt'; 
X_OC(casE) = 0.03;
st{casE} = strcat('TRI - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
casE = casE + 1;

% Vtail
casos{casE} = 'vtail_T1-25_0 m_s-TriUniform-ThickSurf-x0_018m.txt'; 
X_OC(casE) = 0.018;
st{casE} = strcat('Tri - Vee (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
casE = casE + 1;

% Vtail, Sideslip
casos{casE} = 'vtail_T5-a0_0°-25_0m_s-TriUniform-ThickSurf-x0_018m.txt'; 
X_OC(casE) = 0.018;
st{casE} = strcat('Tri - Vee (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
asE = casE + 1;

prefix = strcat('MILVUS_');
mark_legend = st;