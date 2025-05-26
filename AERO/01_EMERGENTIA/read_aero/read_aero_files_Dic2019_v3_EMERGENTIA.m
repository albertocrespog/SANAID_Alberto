function [casos prefix mark_legend X_OC] = read_aero_files_June2019_v2(V_Performance)

V = V_Performance.V;
h = V_Performance.h;

% Estudios
% Front Wing analysis for different locations
casE = 1;
% casos{casE} = 'T1-85_0 m_s-LLT-0_0kg-x0_0m.txt'; % Front Wing alone meassured @ Leading Edge 
% X_OC(casE) = 0.000;
% st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% casE = casE + 1;
% casos{casE} = 'T1-85_0 m_s-LLT-0_0kg-x0_02m.txt'; % Front Wing alone meassured @ -0.005 from Leading Edge
% X_OC(casE) = 0.005;
% st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% casE = casE + 1;
% casos{casE} = 'T1-85_0 m_s-LLT-0_0kg-x0_03m.txt'; % Front Wing alone meassured @ -0.010 from Leading Edge 
% X_OC(casE) = 0.010;
% st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% casE = casE + 1;
% casos{casE} = 'T1-85_0 m_s-LLT-0_0kg-x0_04m.txt'; % Front Wing alone meassured @ -0.015 from Leading Edge 
% X_OC(casE) = 0.015;
% st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% casE = casE + 1;
% casos{casE} = 'T1-85_0 m_s-LLT-0_0kg-x0_05m.txt'; % Front Wing alone meassured @ -0.020 from Leading Edge
% X_OC(casE) = 0.020;
% st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% casE = casE + 1;
% casos{casE} = 'T1-85_0 m_s-LLT-0_0kg-x0_055m.txt'; % Front Wing alone meassured @ -0.025 from Leading Edge 
% X_OC(casE) = 0.025;
% st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% casE = casE + 1;
% casos{casE} = 'T1-85_0 m_s-LLT-0_0kg-x0_06m.txt'; % Front Wing alone meassured @ -0.030 from Leading Edge
% X_OC(casE) = 0.030;
% st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% casE = casE + 1;
% casos{casE} = 'T1-85_0 m_s-LLT-0_0kg-x0_06m.txt'; % Front Wing alone meassured @ -0.030 from Leading Edge
% X_OC(casE) = 0.035;
% st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% casE = casE + 1;
% casos{casE} = 'T1-85_0 m_s-LLT-0_0kg-x0_0m.txt'; % Front Wing alone meassured @ Leading Edge
% X_OC(casE) = 0.040;
% st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% casE = casE + 1;
% casos{casE} = 'T1-85_0 m_s-LLT-0_0kg-x0_05m.txt'; % Front Wing alone meassured @ 0.005 from Leading Edge
% X_OC(casE) = 0.041;
% st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% casE = casE + 1;
% casos{casE} = 'T1-85_0 m_s-LLT-0_0kg-x0_1m.txt'; % Front Wing alone meassured @ 0.010 from Leading Edge
% X_OC(casE) = 0.042;
% st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% casE = casE + 1;
% casos{casE} = 'T1-85_0 m_s-LLT-0_0kg-x0_125m.txt'; % Front Wing alone meassured @ 0.100 from Leading Edge
% X_OC(casE) = 0.043;
% st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% casE = casE + 1;
% casos{casE} = 'T1-85_0 m_s-LLT-0_0kg-x0_13m.txt'; % Front Wing alone meassured @ 0.140 from Leading Edge
% X_OC(casE) = 0.044;
% st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% casE = casE + 1;
% casos{casE} = 'T1-85_0 m_s-LLT-0_0kg-x0_15m.txt'; % Front Wing alone meassured @ 0.150 from Leading Edge
% X_OC(casE) = 0.045;
% st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
% casE = casE + 1;
casos{casE} = 'T1-60_0 m_s-VLM2-0_0kg-x0_338m_W1.txt'; % Front Wing alone meassured @ 0.338 from Leading Edge
X_OC(casE) = 0.338;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
casE = casE + 1;
casos{casE} = 'T1-60_0 m_s-VLM2-0_0kg-x0_129_z_0_15_W2.txt'; % Vtail alone meassured @ 0.129 from Leading Edge
X_OC(casE) = 0.129;
st{casE} = strcat('LLT - W_2 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m');
casE = casE + 1;

prefix = strcat('PepiñoXXXL_');

mark_legend = st;