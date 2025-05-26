function [casos prefix mark_legend X_OC] = read_aero_files_April2024_SOLARTII_v1(Performance)

V = Performance.V;
h = Performance.h;
 
% Wing 11, sweep 14.25, taper = 0.6, dihedral = 0 deg - AC1_Sweep_14_25_diedral_0_ T1-20.0 m/s-LLT-x0.75m-z0.0m_h1000
casE = 1;
% wing
casos{casE} = 'T1-20_0 m_s-TriUniform-ThickSurf-x0_75m-z0_0m_dihedral_0_sweep_14_25.txt'; 
X_OC(casE) = 0.75;
st{casE} = strcat('TriUniform - W_1, \Lambda =',num2str(14.25), '\circ, \lambda=',num2str(0.6),', \Gamma = ',num2str(0),'\circ, X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

% Wing 12, sweep 14.25, taper = 0.6, dihedral = 2.5 deg
casE = casE + 1;
casos{casE} = 'T1-20_0 m_s-TriUniform-ThickSurf-x0_67m-z0_0m_dihedral_2_5_sweep_14_25.txt';
X_OC(casE) = 0.65;
st{casE} = strcat('TriUniform - W_1, \Lambda =',num2str(14.25), '\circ, \lambda=',num2str(0.6),', \Gamma = ',num2str(2.5),'\circ, X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

% Wing 13, sweep 14.25, taper = 0.6, dihedral = 5 deg
casE = casE + 1;
casos{casE} = 'T1-20_0 m_s-TriUniform-ThickSurf-x0_65m-z0_0m_dihedral_5_sweep_14_25.txt';
X_OC(casE) = 0.65;
st{casE} = strcat('TriUniform - W_1, \Lambda =',num2str(14.25), '\circ, \lambda=',num2str(0.6),', \Gamma = ',num2str(5),'\circ, X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

% Wing 21, sweep 25, taper = 0.6, dihedral = 0 deg
casE = casE + 1;
% wing
casos{casE} = 'T1-20_0 m_s-VLM1-x0_82m-z0_0m_dihedral_0_sweep_25.txt'; 
X_OC(casE) = 0.82;
st{casE} = strcat('VLM - W_1, \Lambda =',num2str(25), '\circ, \lambda=',num2str(0.6),', \Gamma = ',num2str(0),'\circ, m X_{ac} =',num2str(X_OC(casE)),' CASE=',num2str(casE));

% Wing 22, sweep 25, taper = 0.6, dihedral = 2.5 deg
casE = casE + 1;
casos{casE} = 'T1-20_0 m_s-TriUniform-ThinSurf-x0_83m-z0_0m_dihedral_2_5_sweep_25.txt';
X_OC(casE) = 0.83;
st{casE} = strcat('TriUniform - W_1, \Lambda =',num2str(25), '\circ, \lambda=',num2str(0.6),', \Gamma = ',num2str(2.5),'\circ, m X_{ac} =',num2str(X_OC(casE)),' CASE=',num2str(casE));

% Wing 23, sweep 25, taper = 0.6, dihedral = 5 deg
casE = casE + 1;
casos{casE} = 'T1-20_0 m_s-TriUniform-ThinSurf-x0_83m-z0_0m_dihedral_5_sweep_25.txt';
X_OC(casE) = 0.83;
st{casE} = strcat('TriUniform - W_1, \Lambda =',num2str(25), '\circ, \lambda=',num2str(0.6),', \Gamma = ',num2str(5),'\circ, m X_{ac} =',num2str(X_OC(casE)),' CASE=',num2str(casE));

% Wing 31, sweep 35, taper = 0.6, dihedral = 0 deg
casE = casE + 1;
% wing
casos{casE} = 'T1-20_0 m_s-VLM1-x1_03m-z0_0m_dihedral_0_sweep_35.txt'; 
X_OC(casE) = 1.03;
st{casE} = strcat('VLM - W_1, \Lambda =',num2str(35), '\circ, \lambda=',num2str(0.6),', \Gamma = ',num2str(0),'\circ, m X_{ac} =',num2str(X_OC(casE)),' CASE=',num2str(casE));

% Wing 32, sweep 35, taper = 0.6, dihedral = 2.5 deg
casE = casE + 1;
casos{casE} = 'T1-20_0 m_s-TriUniform-ThinSurf-x1_03m-z0_0m_dihedral_2_5_sweep_35.txt';
X_OC(casE) = 0.83;
st{casE} = strcat('TriUniform - W_1, \Lambda =',num2str(35), '\circ, \lambda=',num2str(0.6),', \Gamma = ',num2str(2.5),'\circ, m X_{ac} =',num2str(X_OC(casE)),' CASE=',num2str(casE));

% Wing 33, sweep 35, taper = 0.6, dihedral = 5 deg
casE = casE + 1;
casos{casE} = 'T1-20_0 m_s-TriUniform-ThinSurf-x1_03m-z0_0m_dihedral_5_sweep_35.txt';
X_OC(casE) = 0.83;
st{casE} = strcat('TriUniform - W_1, \Lambda =',num2str(35), '\circ, \lambda=',num2str(0.6),', \Gamma = ',num2str(5),'\circ, m X_{ac} =',num2str(X_OC(casE)),' CASE=',num2str(casE));

% Stdudies done initially
casE = casE + 1;
% wing
casos{casE} = 'T1-20_0 m_s-LLT-x0_52m_3000m.txt'; % Case 1
X_OC(casE) = 0.085;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

% VTP Longitudinal
casE = casE + 1;
casos{casE} = 'T5-a0_0-20_0m_s-VLM2-x0_5m.txt'; % Case 2
X_OC(casE) = 0.06475;
st{casE} = strcat('VLM - VTP (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

prefix = strcat('SOLARTII_'); 

mark_legend = st;


