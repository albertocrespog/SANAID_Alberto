function [casos prefix mark_legend X_OC] = read_aero_files_March2022_EMERGENTIA(Performance)

V = Performance.V;
h = Performance.h;

% Estudios
% Front Wing analysis for different locations
casE = 1;
casos{casE} = 'Tandem_wing_T1-26_0 m_s-LLT-x0_0505m.txt'; % Case 1
X_OC(casE) = 0.0505;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
% Canard
casE = casE + 1;
% casos{casE} = 'Tandem_canard_T1-26_0 m_s-LLT-x0_0505m.txt'; % Case 2
casos{casE} = 'WING-T6_1-TriUniform-ThickSurf-x0_05m-DEF.txt'; % Case 2
X_OC(casE) = 0.0505;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
% Vee Tail No center section
casE = casE + 1;
% casos{casE} = 'Vee_projected_NCS_T1-26_0 m_s-LLT-x0_082m.txt';  % Case 3
casos{casE} = 'VTAIL-T6_1-TriUniform-ThickSurf-vl-x0_09m-DEF.txt';  % Case 3
X_OC(casE) = 0.082;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
% Vee Tail center section
casE = casE + 1;
casos{casE} = 'VTAIL-T6_1-TriUniform-ThickSurf-vl-x0_09m-CENTER_SECTION.txt';  % Case 4
X_OC(casE) = 0.071;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
% VTP for tandem configuration
casos{casE} = 'Tandem_W&VTP_T1-26_0 m_s-LLT-x0_0m.txt';  % Case 5
X_OC(casE) = 0.0;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
% CYbeta Vtail file No center section
casos{casE} = 'VTAIL-T6_1-TriUniform-ThickSurf-beta-visc.loop.txt';  % Case 6
X_OC(casE) = 0.0;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
casE = casE + 1;
% CYbeta Vtail file WITH center section
casos{casE} = 'VTAIL-T6_1-TriUniform-ThickSurf-beta-center_section-vl.txt';  % Case 7
X_OC(casE) = 0.0;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

prefix = strcat('EMERGENTIA_'); 

mark_legend = st;