function Weight_mission = T120_weightconfiguration

% excel_file = '../AIRCRAFT/TARSIS/Weights_CG/Tabla_condiciones_CG_T120.xlsx';
excel_file = 'Tabla_condiciones_CG_T120_v4.xlsx';

% Sheet for full fuel 
Sheet1 = cell2mat(readcell(excel_file,'Sheet','CAMARA_E180','Range','D17:H28'));
% Sheet for half fuel 
Sheet2 = cell2mat(readcell(excel_file,'Sheet','CAMARA_E180','Range','I17:M28'));
% Sheet for no fuel 
% Sheet3 = cell2mat(readcell(excel_file,'Sheet','CAMARA_E140','Range','N17:R28'));
Sheet3 = cell2mat(readcell(excel_file,'Sheet','CAMARA_E180','Range','N17:R28'));

%% Data for all fuel 
m_misil = 3;
N_mis_vec_f  = Sheet1(1,:);
M_mis_vec_f  = Sheet1(1,:).*m_misil;
W_vec_f  = Sheet1(2,:);
W_fuel_vec_f = Sheet1(3,:);

% Centros de Gravedad
x_XCG_vec_f  = Sheet1(4,:)./1000;
y_XCG_vec_f  = Sheet1(5,:)./1000;
z_XCG_vec_f  = Sheet1(6,:)./1000;

% inercias
I_xx_vec_f  = Sheet1(7,:);
I_yy_vec_f  = Sheet1(8,:);
I_zz_vec_f  = Sheet1(9,:);

% inercias cruzadas
I_xy_vec_f  = Sheet1(10,:);
I_xz_vec_f  = -Sheet1(11,:);
I_yz_vec_f  = -Sheet1(12,:);

Weight_mission.W_fuel_vec_f = W_fuel_vec_f;
Weight_mission.N_mis_vec_f = N_mis_vec_f;
Weight_mission.M_mis_vec_f = M_mis_vec_f;
Weight_mission.W_vec_f = W_vec_f;

Weight_mission.x_XCG_vec_f = x_XCG_vec_f;
Weight_mission.y_XCG_vec_f = y_XCG_vec_f;
Weight_mission.z_XCG_vec_f = z_XCG_vec_f;

Weight_mission.I_xx_vec_f = I_xx_vec_f;
Weight_mission.I_yy_vec_f = I_yy_vec_f;
Weight_mission.I_zz_vec_f = I_zz_vec_f;

Weight_mission.I_xy_vec_f = I_xy_vec_f;
Weight_mission.I_xz_vec_f = I_xz_vec_f;
Weight_mission.I_yz_vec_f = I_yz_vec_f;

%% Data for half fuel 
m_misil = 3;
N_mis_vec_f2  = Sheet2(1,:);
M_mis_vec_f2  = Sheet2(1,:).*m_misil;
W_fuel_vec_f2 = Sheet2(3,:);
W_vec_f2  = Sheet2(2,:);

% inercias
x_XCG_vec_f2  = Sheet2(4,:)./1000;
y_XCG_vec_f2  = Sheet2(5,:)./1000;
z_XCG_vec_f2  = Sheet2(6,:)./1000;

% inercias
I_xx_vec_f2  = Sheet2(7,:);
I_yy_vec_f2  = Sheet2(8,:);
I_zz_vec_f2  = Sheet2(9,:);

% inercias cruzadas
I_xy_vec_f2  = Sheet2(10,:);
I_xz_vec_f2  = -Sheet2(11,:);
I_yz_vec_f2  = -Sheet2(12,:);

Weight_mission.W_fuel_vec_f2 = W_fuel_vec_f2;
Weight_mission.N_mis_vec_f2 = N_mis_vec_f2;
Weight_mission.M_mis_vec_f2 = M_mis_vec_f2;
Weight_mission.W_vec_f2 = W_vec_f2;

Weight_mission.x_XCG_vec_f2 = x_XCG_vec_f2;
Weight_mission.y_XCG_vec_f2 = y_XCG_vec_f2;
Weight_mission.z_XCG_vec_f2 = z_XCG_vec_f2;

Weight_mission.I_xx_vec_f2 = I_xx_vec_f2;
Weight_mission.I_yy_vec_f2 = I_yy_vec_f2;
Weight_mission.I_zz_vec_f2 = I_zz_vec_f2;

Weight_mission.I_xy_vec_f2 = I_xy_vec_f2;
Weight_mission.I_xz_vec_f2 = I_xz_vec_f2;
Weight_mission.I_yz_vec_f2 = I_yz_vec_f2;

%% Data for NO fuel 
m_misil = 3;
N_mis_vec_nf  = Sheet3(1,:);
M_mis_vec_nf  = Sheet3(1,:).*m_misil;
W_vec_nf  = Sheet3(2,:);
W_fuel_vec_nf = Sheet3(3,:);

%centros de gravedad
x_XCG_vec_nf  = Sheet3(4,:)./1000;
y_XCG_vec_nf  = Sheet3(5,:)./1000;
z_XCG_vec_nf  = Sheet3(6,:)./1000;

% inercias
I_xx_vec_nf  = Sheet3(7,:);
I_yy_vec_nf  = Sheet3(8,:);
I_zz_vec_nf  = Sheet3(9,:);

% inercias cruzadas
I_xy_vec_nf  = Sheet3(10,:);
I_xz_vec_nf  = -Sheet3(11,:);
I_yz_vec_nf  = -Sheet3(12,:);

Weight_mission.W_fuel_vec_nf = W_fuel_vec_nf;
Weight_mission.N_mis_vec_nf = N_mis_vec_nf;
Weight_mission.M_mis_vec_nf = M_mis_vec_nf;
Weight_mission.W_vec_nf = W_vec_nf;

Weight_mission.x_XCG_vec_nf = x_XCG_vec_nf;
Weight_mission.y_XCG_vec_nf = y_XCG_vec_nf;
Weight_mission.z_XCG_vec_nf = z_XCG_vec_nf;

Weight_mission.I_xx_vec_nf = I_xx_vec_nf;
Weight_mission.I_yy_vec_nf = I_yy_vec_nf;
Weight_mission.I_zz_vec_nf = I_zz_vec_nf;

Weight_mission.I_xy_vec_nf = I_xy_vec_nf;
Weight_mission.I_xz_vec_nf = I_xz_vec_nf;
Weight_mission.I_yz_vec_nf = I_yz_vec_nf;
