function [casos prefix mark_legend X_OC] = read_aero_files_June2019(V_Performance)
% Estudios
% Front Wing analysis for different locations
% Front Wing alone meassured @ Leading Edge
% automatic assignment of cases
V = V_Performance.V;
h = V_Performance.h;
% casE = 1;
% casos{1} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_000m.txt'; 
% % casos{casE} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_000m.txt'; 
% X_OC(casE) = 0.000;
% st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% casE = casE + 1;
% casos{casE} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_005m.txt'; 
% X_OC(casE) = 0.005;
% st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% casE = casE + 1;
% casos{casE} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_010m.txt'; 
% X_OC(casE) = 0.010;
% st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% casE = casE + 1;
% casos{casE} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_015m.txt'; 
% X_OC(casE) = 0.015;
% st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% casE = casE + 1;
% casos{casE} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_020m.txt'; 
% X_OC(casE) = 0.020;
% st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% casE = casE + 1;
% casos{casE} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_025m.txt'; 
% X_OC(casE) = 0.025;
% st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% casE = casE + 1;
% casos{casE} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_030m.txt'; 
% X_OC(casE) = 0.030;
% st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% casE = casE + 1;
% casos{casE} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_035m.txt'; 
% X_OC(casE) = 0.035;
% st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% casE = casE + 1;
% casos{casE} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_040m.txt'; 
% X_OC(casE) = 0.040;
% st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% casE = casE + 1;
% casos{casE} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_044m.txt'; 
% X_OC(casE) = 0.044;
% st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% casE = casE + 1;
% casos{casE} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_045m.txt'; 
% X_OC(casE) = 0.045;
% st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% casE = casE + 1;
% casos{casE} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_046m.txt'; 
% X_OC(casE) = 0.046;
% st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% casE = casE + 1;
casos{1} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_000m.txt'; % Front Wing alone meassured @ -0.005 from Leading Edge
casos{2} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_005m.txt'; % Front Wing alone meassured @ -0.005 from Leading Edge
casos{3} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_010m.txt'; % Front Wing alone meassured @ -0.010 from Leading Edge
casos{4} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_015m.txt'; % Front Wing alone meassured @ -0.015 from Leading Edge
casos{5} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_020m.txt'; % Front Wing alone meassured @ -0.020 from Leading Edge
casos{6} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_025m.txt'; % Front Wing alone meassured @ -0.025 from Leading Edge
casos{7} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_030m.txt'; % Front Wing alone meassured @ -0.030 from Leading Edge
casos{8} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_035m.txt'; % Front Wing alone meassured @ -0.031 from Leading Edge
casos{9} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_040m.txt'; % Front Wing alone meassured @ -0.031 from Leading Edge
casos{10} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_044m.txt'; % Front Wing alone meassured @ -0.031 from Leading Edge
casos{11} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_045m.txt'; % Front Wing alone meassured @ -0.031 from Leading Edge
casos{12} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_046m.txt'; % Front Wing alone meassured @ -0.031 from Leading Edge

% Value at which have been meassured for calculating moments about a
% different point

X_OC(1) = 0.000;
X_OC(2) = 0.005;
X_OC(3) = 0.010;
X_OC(4) = 0.015;
X_OC(5) = 0.020;
X_OC(6) = 0.025;
X_OC(7) = 0.030;
X_OC(8) = 0.035;
X_OC(9) = 0.040;
X_OC(10) = 0.044;
X_OC(11) = 0.045;
X_OC(12) = 0.046;

% % Rear WingFront Wing (V-Tail) analysis for different locations
% casos{13} = 'VTOL_SINGLE_VTail_0008_LAST_T1-25_0 m_s-LLT-0_0kg-x0_000m.txt'; % Front Wing alone meassured @ Leading Edge
% casos{14} = 'VTOL_SINGLE_VTail_0008_LAST_T1-25_0 m_s-LLT-0_0kg-x0_005m.txt'; % Front Wing alone meassured @ 0.005 from Leading Edge
% casos{15} = 'VTOL_SINGLE_VTail_0008_LAST_T1-25_0 m_s-LLT-0_0kg-x0_010m.txt'; % Front Wing alone meassured @ 0.010 from Leading Edge
% casos{16} = 'VTOL_SINGLE_VTail_0008_LAST_T1-25_0 m_s-LLT-0_0kg-x0_100m.txt'; % Front Wing alone meassured @ 0.100 from Leading Edge
% casos{17} = 'VTOL_SINGLE_VTail_0008_LAST_T1-25_0 m_s-LLT-0_0kg-x0_140m.txt'; % Front Wing alone meassured @ 0.140 from Leading Edge
% casos{18} = 'VTOL_SINGLE_VTail_0008_LAST_T1-25_0 m_s-LLT-0_0kg-x0_150m.txt'; % Front Wing alone meassured @ 0.150 from Leading Edge

prefix = strcat('ProVANT4.0_');

V = V_Performance.V;
h = V_Performance.h;
% Generates the Legend Dymensions according tot he number of test sts
% Front Wing Legends
st1 = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} = 0,000m');
st2 = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m, X_{ac} = 0.005m');
st3 = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m, X_{ac} = 0.010m');
st4 = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m, X_{ac} = 0.015m');
st5 = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m, X_{ac} = 0.020m');
st6 = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m, X_{ac} = 0.025m');
st7 = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m, X_{ac} = 0.030m');
st8 = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m, X_{ac} = 0.035m');
st9 = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m, X_{ac} = 0.040m');
st10 = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m, X_{ac} = 0.044m');
st11 = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m, X_{ac} = 0.045m');
st12 = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m, X_{ac} = 0.046m');

% % Rear Wing Legends
% st9 = strcat('LLT - V_{tail} (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} = 0,000 m');
% st10 = strcat('LLT - V_{tail} (',num2str(V), 'm/s) @',num2str(h),' m, X_{ac} = 0.005m');
% st11 = strcat('LLT - V_{tail} (',num2str(V), 'm/s) @',num2str(h),' m, X_{ac} = 0.010m');
% st12 = strcat('LLT - V_{tail} (',num2str(V), 'm/s) @',num2str(h),' m, X_{ac} = 0.100m');
% st13 = strcat('LLT - V_{tail} (',num2str(V), 'm/s) @',num2str(h),' m, X_{ac} = 0.140m');
% st14 = strcat('LLT - V_{tail} (',num2str(V), 'm/s) @',num2str(h),' m, X_{ac} = 0.150m');

% Creates the vector with all the legends
% mark_legend = {st1,st2,st3,st4,st5,st6,st7,st8,st9,st10,st11,st12,st13,st14};
mark_legend = {st1,st2,st3,st4,st5,st6,st7,st8,st9,st10,st11,st12};