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
casos{casE} = 'T6_1-TriUniform-ThickSurf-alpha-dihedral30.txt';  % Case 4
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
casos{casE} = 'T6_1-TriUniform-ThickSurf-beta-dihedral30.txt';  % Case 7
X_OC(casE) = 0.0;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

%%w2_dihedral_variation
% Vee Tail center section
casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-alpha-dihedral35.txt';  % Case 8
X_OC(casE) = 0.071;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
% Vee Tail center section
casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-alpha-dihedral40.txt';  % Case 9
X_OC(casE) = 0.071;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
% Vee Tail center section
casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-alpha-dihedral45.txt';  % Case 10
X_OC(casE) = 0.071;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
% Vee Tail center section
casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-alpha-dihedral50.txt';  % Case 11
X_OC(casE) = 0.071;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
% Vee Tail center section
casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-alpha-dihedral55.txt';  % Case 12
X_OC(casE) = 0.071;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));


% CYbeta Vtail file WITH center section
casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-beta-dihedral35.txt';  % Case 13
X_OC(casE) = 0.0;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
% CYbeta Vtail file WITH center section
casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-beta-dihedral40.txt';  % Case 14
X_OC(casE) = 0.0;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
% CYbeta Vtail file WITH center section
casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-beta-dihedral45.txt';  % Case 15
X_OC(casE) = 0.0;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
% CYbeta Vtail file WITH center section
casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-beta-dihedral50.txt';  % Case 16
X_OC(casE) = 0.0;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
% CYbeta Vtail file WITH center section
casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-beta-dihedral55.txt';  % Case 17
X_OC(casE) = 0.0;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

% Vee Tail center section 45º VEE DIHEDRAL, S_h_proj = 30º vee dihedral
casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-alpha-dihedral45_S_Proj.txt';  % Case 18
X_OC(casE) = 0.071;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));
% CYbeta Vtail file WITH center section
casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-beta-dihedral45_S_Proj.txt';  % Case 19
X_OC(casE) = 0.0;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));


% Vtail geometry variation. Cases 1, 2a & 2b files:

%% alpha_var

%dihedral 35 degrees
casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-alpha-x0_0715m-35-case1.txt';  % Case 20
X_OC(casE) = 0.0715;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-alpha-x0_063m-35-case2a.txt';  % Case 21
X_OC(casE) = 0.063;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-alpha-x0_083m-35-case2b.txt';  % Case 22
X_OC(casE) = 0.083;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

%dihedral 40 degrees

casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-alpha-x0_07m-40-case1.txt';  % Case 23
X_OC(casE) = 0.07;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-alpha-x0_067m-40-case2a.txt';  % Case 24
X_OC(casE) = 0.067;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-alpha-x0_075m-40-case2b.txt';  % Case 25
X_OC(casE) = 0.075;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

%dihedral 45 degrees

casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-alpha-x0_065-45-nominal.txt';  % Case 26
X_OC(casE) = 0.065;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

%dihedral 50 degrees

casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-alpha-x0_064m-50-case1.txt';  % Case 27
X_OC(casE) = 0.064;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-alpha-x0_068m-50-case2a.txt';  % Case 28
X_OC(casE) = 0.068;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-alpha-x0_06m-50-case2b.txt';  % Case 29
X_OC(casE) = 0.060;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

%dihedral 55 degrees

casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-alpha-x0_06m-55-case1.txt';  % Case 30
X_OC(casE) = 0.060;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-alpha-x0_069m-55-case2a.txt';  % Case 31
X_OC(casE) = 0.069;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-alpha-x0_057m-55-case2b.txt';  % Case 32
X_OC(casE) = 0.057;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));


%% beta_var

%dihedral 35 degrees
casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-beta-35-case1.txt';  % Case 33
X_OC(casE) = 0.0715;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-beta-35-case2a.txt';  % Case 34
X_OC(casE) = 0.0715;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-beta-35-case2b.txt';  % Case 35
X_OC(casE) = 0.0715;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

%dihedral 40 degrees

casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-beta-40-case1.txt';  % Case 36
X_OC(casE) = 0.0715;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-beta-40-case2a.txt';  % Case 37
X_OC(casE) = 0.0715;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-beta-40-case2b.txt';  % Case 38
X_OC(casE) = 0.0715;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

%dihedral 45 degrees

casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-beta-45-nominal.txt';  % Case 39
X_OC(casE) = 0.0715;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));



%dihedral 50 degrees

casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-beta-50-case1.txt';  % Case 40
X_OC(casE) = 0.0715;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-beta-50-case2a.txt';  % Case 41
X_OC(casE) = 0.0715;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-beta-50-case2b.txt';  % Case 42
X_OC(casE) = 0.0715;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

%dihedral 55 degrees
casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-beta-55-case1.txt';  % Case 43
X_OC(casE) = 0.0715;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-beta-55-case2a.txt';  % Case 44
X_OC(casE) = 0.0715;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));

casE = casE + 1;
casos{casE} = 'T6_1-TriUniform-ThickSurf-beta-55-case2b.txt';  % Case 45
X_OC(casE) = 0.0715;
st{casE} = strcat('LLT - W_1 (',num2str(V), 'm/s) @',num2str(h),' m X_{ac} =',num2str(X_OC(casE)),' m,CASE=',num2str(casE));


















prefix = strcat('EMERGENTIA_'); 

mark_legend = st;