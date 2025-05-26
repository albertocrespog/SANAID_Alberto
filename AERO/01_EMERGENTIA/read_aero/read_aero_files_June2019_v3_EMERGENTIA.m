function [casos prefix mark_legend X_OC] = read_aero_files_June2019_v3_EMERGENTIA(Performance)

V = Performance.V;
h = Performance.h;

% Estudios
% Front Wing analysis for different locations
casE = 1;
casos{casE} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_000m.txt'; 
X_OC(casE) = 0.000;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
casos{casE} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_005m.txt'; 
X_OC(casE) = 0.005;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
casos{casE} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_010m.txt'; 
X_OC(casE) = 0.010;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
casos{casE} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_015m.txt'; 
X_OC(casE) = 0.015;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
casos{casE} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_020m.txt'; 
X_OC(casE) = 0.020;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
casos{casE} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_025m.txt'; 
X_OC(casE) = 0.025;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
casos{casE} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_030m.txt'; 
X_OC(casE) = 0.030;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
casos{casE} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_035m.txt'; 
X_OC(casE) = 0.035;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
casos{casE} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_040m.txt'; 
X_OC(casE) = 0.040;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
casos{casE} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_041m.txt'; 
X_OC(casE) = 0.041;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
casos{casE} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_042m.txt'; 
X_OC(casE) = 0.042;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
casos{casE} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_043m.txt'; 
X_OC(casE) = 0.043;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
casos{casE} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_044m.txt'; 
X_OC(casE) = 0.044;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
casos{casE} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_045m.txt'; 
X_OC(casE) = 0.045;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
casos{casE} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_046m.txt'; 
X_OC(casE) = 0.046;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
casos{casE} = 'EMERGENTIA_Sw1_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_047m.txt'; 
X_OC(casE) = 0.047;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;

% Rear WingFront Wing (V-Tail) analysis for different locations
casos{casE} = 'EMERGENTIA_Sw2_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_000m.txt'; 
X_OC(casE) = 0.000;
st{casE} = strcat('LLT - W_2 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
casos{casE} = 'EMERGENTIA_Sw2_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_005m.txt'; 
X_OC(casE) = 0.005;
st{casE} = strcat('LLT - W_2 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
casos{casE} = 'EMERGENTIA_Sw2_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_010m.txt'; 
X_OC(casE) = 0.010;
st{casE} = strcat('LLT - W_2 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
casos{casE} = 'EMERGENTIA_Sw2_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_015m.txt'; 
X_OC(casE) = 0.015;
st{casE} = strcat('LLT - W_2 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
casos{casE} = 'EMERGENTIA_Sw2_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_020m.txt'; 
X_OC(casE) = 0.020;
st{casE} = strcat('LLT - W_2 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
casos{casE} = 'EMERGENTIA_Sw2_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_025m.txt'; 
X_OC(casE) = 0.025;
st{casE} = strcat('LLT - W_2 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
casos{casE} = 'EMERGENTIA_Sw2_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_030m.txt'; 
X_OC(casE) = 0.030;
st{casE} = strcat('LLT - W_2 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
casos{casE} = 'EMERGENTIA_Sw2_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_035m.txt'; 
X_OC(casE) = 0.035;
st{casE} = strcat('LLT - W_2 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
casos{casE} = 'EMERGENTIA_Sw2_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_036m.txt'; 
X_OC(casE) = 0.036;
st{casE} = strcat('LLT - W_2 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
casos{casE} = 'EMERGENTIA_Sw2_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_037m.txt'; 
X_OC(casE) = 0.037;
st{casE} = strcat('LLT - W_2 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
casos{casE} = 'EMERGENTIA_Sw2_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_038m.txt'; 
X_OC(casE) = 0.038;
st{casE} = strcat('LLT - W_2 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
casos{casE} = 'EMERGENTIA_Sw2_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_039m.txt'; 
X_OC(casE) = 0.039;
st{casE} = strcat('LLT - W_2 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
casos{casE} = 'EMERGENTIA_Sw2_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_040m.txt'; 
X_OC(casE) = 0.040;
st{casE} = strcat('LLT - W_2 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
casos{casE} = 'EMERGENTIA_Sw2_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_041m.txt'; 
X_OC(casE) = 0.041;
st{casE} = strcat('LLT - W_2 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
casos{casE} = 'EMERGENTIA_Sw2_e_dihedral_T1-25_0 m_s-LLT-proj_area_x0_045m.txt'; 
X_OC(casE) = 0.045;
st{casE} = strcat('LLT - W_2 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;

% Modificado por Kiko Almagro. Implementación de la geometría proyectada de la cola en V.
% Aplicación de la teoría del NACA Report. No. 823.
% 4 casos obtenidos en XFLR5: LLT, VLM1, VLM2, Paneles. Se recomienda el
% uso del método VLM2. Sin cajón central.
% 09/11/2021.
% V-tail
casos{casE} = 'T1-30_0 m_s-LLT_vtail_projeted_no_center.txt'; %index_w2 = 32;
X_OC(casE) = 0.030;
st{casE} = strcat('LLT - W_2 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;

casos{casE} = 'T1-30_0 m_s-VLM1_vtail_projected_no_center.txt';%index_w2 = 33;
X_OC(casE) = 0.030;
st{casE} = strcat('VLM1 - W_2 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;

casos{casE} = 'T1-30_0 m_s-VLM2_vtail_projected_no_center.txt';  %index_w2 = 34;
X_OC(casE) = 0.030;
st{casE} = strcat('VLM2 - W_2 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;

casos{casE} = 'T1-30_0 m_s-Panel_vtail_projected_no_center.txt'; %index_w2 = 35;
X_OC(casE) = 0.030;
st{casE} = strcat('Panels - W_2 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1; 
% wing
casos{casE} = 'T1-30_0 m_s-VLM2_wing_dihedral.txt'; %index_w2 = 36;
X_OC(casE) = 0.030;
st{casE} = strcat('Panels - W_2 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
casos{casE} = 'T1-30_0 m_s-VLM2-0_0kg-x0_06m-TG.txt'; %index_w2 = 37;
X_OC(casE) = 0.030;
st{casE} = strcat('Panels - W_2 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;

prefix = strcat('EMERGENTIA_'); 

% VTP
casos{casE} = 'HTP_VTP_T5-a0_0°-236_5m_s-Quads-ThinSurfHTPVTP.txt';  %index_w3 = 38;
X_OC(casE) = 0.00;
st{casE} = strcat('LLT - VTP (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;


mark_legend = st;