function [casos prefix mark_legend X_OC] = read_aero_files_March2022_TARSIS(Performance)

V = Performance.V;
h = Performance.h;
 
% Wing Complete alone
casE = 1;
% casos{casE} = 'wing_e_T75_T1-43_2 m_s-LLT-x0_085m.txt'; % Case 1 exposed
% wing
casos{casE} = 'W1_T1-46_9 m_s-LLT-x0_085m.txt'; % Case 1
X_OC(casE) = 0.085;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

% HTP & VTP Longitudinal
casE = casE + 1;
casos{casE} = 'HTP&VTP_T75_T1-29_5 m_s-LLT-x1_5m.txt'; % Case 2
casos{casE} = 'HTP&VTP_T1-46_9 m_s-LLT-x0_155m.txt'; % Case 2
X_OC(casE) = 0.06475;
st{casE} = strcat('LLT - HTP & VTP (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

% HTP & VTP Lateral-directional
casE = casE + 1;
casos{casE} = 'HTP&VTP_T75_T5-a0_0-29_5m_s-VLM2-x0_15m.txt';  % Case 3
X_OC(casE) = 0.06475;
st{casE} = strcat('LLT - HTP & VTP (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
% Vee Tail center section

% % Wing + 20 cm span on each side
% casE = casE + 1;
% casos{casE} = 'W1_T1-46_9 m_s-LLT-x0_085m_b_plusDelta20cm.txt'; % Case 1
% X_OC(casE) = 0.085;
% 
% % Wing + 20 cm span on each side
% casE = casE + 1;
% casos{casE} = 'W1_T1-46_9 m_s-LLT-x0_085m_b_plusDelta30cm.txt'; % Case 1
% X_OC(casE) = 0.085;
% 
% % Wing + 20 cm span on each side
% casE = casE + 1;
% casos{casE} = 'W1_T1-46_9 m_s-LLT-x0_085m_b_plusDelta50cm.txt'; % Case 1
% X_OC(casE) = 0.085;

% Modified Wings 8/4/2022
% T95 wing + 20 cm on each half
% ID_01. c_nom, b+20	
casE = casE + 1;
casos{casE} = 'ID_01_T1-46_9 m_s-LLT-x0_085m_b_plusDelta20cm_wing_c_nom.txt'; 
X_OC(casE) = 0.085;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
% casE = casE + 1;

% T95 wing + 30 cm on each half
% ID_02. c_nom, b+30	
casE = casE + 1;
casos{casE} = 'ID_02_T1-46_9 m_s-LLT-x0_085m_b_plusDelta30cm_wing_c_nom.txt'; 
X_OC(casE) = 0.085;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

% T95 wing + 50 cm on each half
% ID_03. c_nom, b+50	
casE = casE + 1;
casos{casE} = 'ID_03_T1-46_9 m_s-LLT-x0_085m_b_plusDelta50cm_wing_c_nom.txt'; 
X_OC(casE) = 0.085;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

% ID_04. c_a, b+30	
casE = casE + 1;
casos{casE} = 'ID_04_T1-46_9 m_s-LLT-x0_085m_b_plusDelta30cm_wing_ca.txt'; 
X_OC(casE) = 0.085;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

% ID_05. c_a, b+50	
casE = casE + 1;
casos{casE} = 'ID_05_T1-46_9 m_s-LLT-x0_085m_b_plusDelta50cm_wing_ca.txt'; 
X_OC(casE) = 0.085;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

% ID_06. c_a, b+70	
casE = casE + 1;
casos{casE} = 'ID_06_T1-46_9 m_s-LLT-x0_085m_b_plusDelta70cm_wing_ca.txt'; 
X_OC(casE) = 0.085;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
 
% ID_07. c_a, b+50	
casE = casE + 1;
casos{casE} = 'ID_07_T1-46_9 m_s-LLT-x0_085m_b_plusDelta50cm_wing_cb.txt'; 
X_OC(casE) = 0.085;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

% ID_08. c_a, b+70	
casE = casE + 1;
casos{casE} = 'ID_08_T1-46_9 m_s-LLT-x0_085m_b_plusDelta70cm_wing_cb.txt'; 
X_OC(casE) = 0.085;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
 
% ID_09. c_a, b+90	
casE = casE + 1;
casos{casE} = 'ID_09_T1-46_9 m_s-LLT-x0_085m_b_plusDelta90cm_wing_cb.txt'; 
X_OC(casE) = 0.085;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

prefix = strcat('TARSIS_'); 

mark_legend = st;