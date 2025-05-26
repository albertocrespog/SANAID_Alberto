function [casos prefix mark_legend X_OC] = read_aero_files_June2019_v3_PEPINOXXL(V_Performance)

V = V_Performance.V;
h = V_Performance.h;

% Estudios
% Front Wing analysis for different locations
casE = 1;
%casos{casE} = 'PEPINOXXL_W1_T1-40_0 m_s-LLT-x0_325m_flaperon.txt'; 
casos{casE} = 'PEPINOXXL_w1_T1-40_0 m_s-LLT-x0_325m.txt';
X_OC(casE) = 0.325;
st{casE} = strcat('VLM2- W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% Vtail analysis for different locations
casE = casE + 1;
casos{casE} = 'PEPINOXXL_Vee_T1-50_0 m_s-LLT-x0_117m-z0_2m.txt'; % Vtail alone meassured @ 0.129 from Leading Edge
X_OC(casE) = 0.117;
st{casE} = strcat('LLT - W_2 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% Front wing with flaperon
casE = casE + 1;
casos{casE} = 'PEPINOXXL_W1_T1-40_0 m_s-LLT-x0_325m_flaperon.txt'; % Flaperon
X_OC(casE) = 0.325;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% Front wing with flaperon
casE = casE + 1;
casos{casE} = 'PEPINOXXL_W1_T1-40_0 m_s-LLT-x0_325m_flap.txt'; % Flap
X_OC(casE) = 0.325;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% Front wing with flaperon
casE = casE + 1;
casos{casE} = 'PEPINOXXL_W1_T1-40_0 m_s-LLT-x0_325m_flap_flaperon.txt'; % Flaperon y flaps
X_OC(casE) = 0.325;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
prefix = strcat('PepiñoXXXL_');
mark_legend = st;