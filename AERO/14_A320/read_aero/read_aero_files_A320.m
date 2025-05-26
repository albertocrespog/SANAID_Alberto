function [casos prefix mark_legend X_OC] = read_aero_files_A320(Performance)

V = Performance.V;
h = Performance.h;

% Estudios
% Front Wing analysis for different locations
casE = 1;
casos{casE} = 'Wing_M_0_2_T1-66_5 m_s-LLT-x3_875m.txt'; 
casos{casE} = 'WING-T6_1-TriUniform-ThickSurf-x3.85-visc_loop.txt';
X_OC(casE) = 3.875;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;

% HTP
casos{casE} = 'HTPVTP_M_0_2_T1-66_5 m_s-LLT-x4_55m.txt'; 
casos{casE} = 'HTP-T6_1-TriUniform-ThickSurf-x1.7m-case1.txt';
X_OC(casE) = 4.55;
st{casE} = strcat('LLT - HTP (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;

% VTP
casos{casE} = 'VTP_M_0_2_T5-a0_0deg-66_5m_s-VLM2-x4_5m.txt'; 
casos{casE} = 'T6_1-TriUniform-ThickSurf-beta-VTP_A320.txt'; 
X_OC(casE) = 0.00;
st{casE} = strcat('LLT - VTP (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;

% HTP Modifications LLT
% Case 1
casos{casE} = 'T1-99_8 m_s-LLT-x2_11m-case1.txt'; 
casos{casE} = 'HTP-T6_1-TriUniform-ThickSurf-x1.7m-case1.txt';
X_OC(casE) = 2.11;
st{casE} = strcat('LLT - HTP (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;

% Case 2
casos{casE} = 'T1-99_8 m_s-LLT-x2_06m-case2.txt'; 
casos{casE} = 'HTP-T6_1-TriUniform-ThickSurf-x1_6m-visc_loop-case2.txt';
X_OC(casE) = 2.06;
st{casE} = strcat('LLT - HTP (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;

% Case 3
casos{casE} = 'T1-99_8 m_s-LLT-x2_1695m-case3.txt'; 
casos{casE} = 'HTP-T6_1-TriUniform-ThickSurf-x1_8m-visc_loop-case3.txt';

X_OC(casE) = 2.1695;
st{casE} = strcat('LLT - HTP (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;

% Case 4
casos{casE} = 'T1-99_8 m_s-LLT-x0_748m-case4.txt'; 
casos{casE} = 'HTP-T6_1-TriUniform-ThickSurf-x0_65m-visc_loop-case4.txt';
X_OC(casE) = 0.748;
st{casE} = strcat('LLT - HTP (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;

% Case 5
casos{casE} = 'T1-99_8 m_s-LLT-x0_835m-case5.txt'; 
casos{casE} = 'HTP-T6_1-TriUniform-ThickSurf-x0_72m_visc_loop-case5.txt';
X_OC(casE) = 0.835;
st{casE} = strcat('LLT - HTP (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;

% Case 6
casos{casE} = 'T1-99_8 m_s-LLT-x0_6825m-case6.txt'; 
casos{casE} = 'HTP-T6_1-TriUniform-ThickSurf-x0_6m_visc_loop_case6.txt';
X_OC(casE) = 0.6825;
st{casE} = strcat('LLT - HTP (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;

% Case 7
casos{casE} = 'T1-99_8 m_s-LLT-xm0_6145m-case7.txt'; 
casos{casE} = 'HTP-T6_1-TriUniform-ThickSurf-xm0_2m_visc_loop-case7.txt';
X_OC(casE) = -0.6145;
st{casE} = strcat('LLT - HTP (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;

% Case 8
casos{casE} = 'T1-99_8 m_s-LLT-xm0_385m-case8.txt'; 
casos{casE} = 'HTP-T6_1-TriUniform-ThickSurf-x0m-visc_loop-case8.txt';
X_OC(casE) = -0.385;
st{casE} = strcat('LLT - HTP (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;

% Case 9
casos{casE} = 'T1-99_8 m_s-LLT-xm0_805m-case9.txt'; 
casos{casE} = 'HTP-T6_1-TriUniform-ThickSurf-xm0_4_visc_loop-case9.txt';
X_OC(casE) = -0.805;
st{casE} = strcat('LLT - HTP (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;

%% VLM Studies
% Front Wing analysis for different locations
casos{casE} = 'W1_T1-99_8 m_s-VLM2-x4_73m.txt'; 
X_OC(casE) = 4.73;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;

% HTP Modifications LLT
% Case 1
casos{casE} = 'HTP_T1-99_8 m_s-VLM2-x2_11m-case1.txt'; 
X_OC(casE) = 2.11;
st{casE} = strcat('LLT - HTP (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;

% Case 2
casos{casE} = 'HTP_T1-99_8 m_s-VLM2-x2_06m-case2.txt'; 
X_OC(casE) = 2.06;
st{casE} = strcat('LLT - HTP (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;

% Case 3
casos{casE} = 'HTP_T1-99_8 m_s-VLM2-x2_1695m-case3.txt'; 
X_OC(casE) = 2.1695;
st{casE} = strcat('LLT - HTP (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;

% Case 4
casos{casE} = 'HTP_T1-99_8 m_s-VLM2-x0_748m-case4.txt'; 
X_OC(casE) = 0.748;
st{casE} = strcat('LLT - HTP (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;

% Case 5
casos{casE} = 'HTP_T1-99_8 m_s-VLM2-x0_835m-case5.txt'; 
X_OC(casE) = 0.835;
st{casE} = strcat('LLT - HTP (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;

% Case 6
casos{casE} = 'HTP_T1-99_8 m_s-VLM2-x0_6825m-case6.txt'; 
X_OC(casE) = 0.6825;
st{casE} = strcat('LLT - HTP (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;

% Case 7
casos{casE} = 'HTP_T1-99_8 m_s-VLM2-xm0_6145m-case7.txt'; 
X_OC(casE) = -0.6145;
st{casE} = strcat('LLT - HTP (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;

% Case 8
casos{casE} = 'HTP_T1-99_8 m_s-VLM2-xm0_385m-case8.txt'; 
X_OC(casE) = -0.385;
st{casE} = strcat('LLT - HTP (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;

% Case 9
casos{casE} = 'HTP_T1-99_8 m_s-VLM2-xm0_805m-case9.txt'; 
X_OC(casE) = -0.805;
st{casE} = strcat('LLT - HTP (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;

prefix = strcat('A320_'); 

mark_legend = st;