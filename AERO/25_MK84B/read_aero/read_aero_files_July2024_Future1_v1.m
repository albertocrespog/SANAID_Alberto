function [casos prefix mark_legend X_OC] = read_aero_files_MK84_v1(Performance)

V = V_Performance.V;
h = V_Performance.h;

% Estudios
% 1 - Front Wing analysis for different locations
casE = 1;
%casos{casE} = 'PEPINOXXL_W1_T1-40_0 m_s-LLT-x0_325m_flaperon.txt'; 
casos{casE} = 'PEPINOXXL_w1_T1-40_0 m_s-LLT-x0_325m.txt';
X_OC(casE) = 0.325;
st{casE} = strcat('VLM2- W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% 2 - Vtail NACA 0007 analysis for different locations
casE = casE + 1;
casos{casE} = 'PEPINOXXL_Vee_T1-50_0 m_s-VLM1-x0_117m-z0_2m_naca0007.txt'; % Vtail alone meassured @ 0.129 from Leading Edge
X_OC(casE) = 0.117;
st{casE} = strcat('LLT - W_2 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% 3 - Front wing with flaperon
casE = casE + 1;
casos{casE} = 'PEPINOXXL_W1_T1-40_0 m_s-LLT-x0_325m_flaperon.txt'; % Flaperon
X_OC(casE) = 0.325;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% 4 - Front wing with flap
casE = casE + 1;
casos{casE} = 'PEPINOXXL_W1_T1-40_0 m_s-LLT-x0_325m_flap.txt'; % Flap
X_OC(casE) = 0.325;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% 5 - Front wing with flaperon + flap
casE = casE + 1;
casos{casE} = 'PEPINOXXL_W1_T1-40_0 m_s-LLT-x0_325m_flap_flaperon.txt'; % Flaperon y flaps
X_OC(casE) = 0.325;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% 6 - Vtail NACA 0007 analysis for different locations VLM
casE = casE + 1;
casos{casE} = 'PEPINOXXL_Vee_T1-50_0 m_s-VLM1-x0_117m-z0_2m_naca0007.txt'; % Vtail alone meassured @ 0.129 from Leading Edge
X_OC(casE) = 0.117;
% 7 - Vtail NACA 0014 analysis for different locations VLM
casE = casE + 1;
casos{casE} = 'PEPINOXXL_Vee_T1-50_0 m_s-VLM1-x0_117m-z0_2m_naca0014.txt'; % Vtail alone meassured @ 0.129 from Leading Edge
X_OC(casE) = 0.117;
% 8 - Front wing with STCYR
casE = casE + 1;
casos{casE} = 'PEPINOXXL_w1_T1-40_0 m_s-LLT-x0_325m_STCYR.txt'; % Vtail alone meassured @ 0.129 from Leading Edge
X_OC(casE) = 0.325;
% 9 - Front wing with STCYR
casE = casE + 1;
casos{casE} = 'PEPINOXXL_w1_T1-40_0 m_s-LLT-x0_325m_MH32.txt'; % Vtail alone meassured @ 0.129 from Leading Edge
X_OC(casE) = 0.325;
% 10 - Front wing with TIP	
casE = casE + 1;
casos{casE} = 'PEPINOXXL_w1_T1-40_0 m_s-LLT-x0_325m_tip.txt'; % Vtail alone meassured @ 0.129 from Leading Edge
X_OC(casE) = 0.325;
% 11 - Front wing with WINGTIP	
casE = casE + 1;
casos{casE} = 'PEPINOXXL_w1_T1-40_0 m_s-LLT-x0_325m_wingtip.txt'; % Vtail alone meassured @ 0.129 from Leading Edge
X_OC(casE) = 0.325;
% 12 - VTAIL new  NACA0014 	
casE = casE + 1;
casos{casE} = 'PEPINOXXL_Vee_T1-50_0 m_s-VLM1-x0_11m-z0_13m_naca007.txt'; % Vtail alone meassured @ 0.129 from Leading Edge
X_OC(casE) = 0.11;

prefix = strcat('PepiñoXXXL_');
mark_legend = st;