function [casos prefix mark_legend X_OC] = read_aero_files_MILVUS_v3(Performance)

V = Performance.V;
h = Performance.h;

% Estudios
% Front Wing analysis for different locations
casE = 1;
casos{casE} = 'wing_T1-25_0 m_s-TriUniform-ThickSurf-x0_031m.txt'; 
X_OC(casE) = 0.031;
st{casE} = strcat('TRI - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
casE = casE + 1;

% Front Wing analysis SideSlip (same as front)
casos{casE} = 'wing_T5-a0_0°-25_0m_s-TriUniform-ThickSurf-x0_031m.txt'; 
X_OC(casE) = 0.03;
st{casE} = strcat('TRI - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
casE = casE + 1;

% Vtail
casos{casE} = 'vtail_T1-25_0 m_s-TriUniform-ThickSurf-x0_016m.txt'; 
X_OC(casE) = 0.016;
st{casE} = strcat('Tri - Vee (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
casE = casE + 1;

% Vtail, Sideslip
casos{casE} = 'vtail_T5-a0_0°-25_0m_s-TriUniform-ThickSurf-x0_016m.txt'; 
X_OC(casE) = 0.016;
st{casE} = strcat('Tri - Vee (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
casE = casE + 1;

% Extended Wing analysis for different locations
casos{casE} = 'wing_extended_T1-25_0 m_s-TriUniform-ThickSurf-x0_029m.txt'; 
X_OC(casE) = 0.029;
st{casE} = strcat('TRI - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
casE = casE + 1;

% Front Wing analysis SideSlip (same as front)
casos{casE} = 'wing_extended_T5-a0_0°-25_0m_s-TriUniform-ThickSurf-x0_029m.txt'; 
X_OC(casE) = 0.029;
st{casE} = strcat('TRI - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
casE = casE + 1;

% Extended Vtail analysis for different locations
casos{casE} = 'vtail_extended_T1-25_0 m_s-TriUniform-ThickSurf-x0_016m.txt'; 
X_OC(casE) = 0.016;
st{casE} = strcat('TRI - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
casE = casE + 1;

% Front Wing analysis SideSlip (same as front)
casos{casE} = 'vtail_extended_T5-a0_0°-25_0m_s-TriUniform-ThickSurf-x0_0m.txt'; 
X_OC(casE) = 0.016;
st{casE} = strcat('TRI - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
casE = casE + 1;



prefix = strcat('MILVUS_');
mark_legend = st;