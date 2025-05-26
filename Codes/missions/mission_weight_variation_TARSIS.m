function [conditions OUTPUT_read_XLSX] = mission_weight_variation_TARSIS(OUTPUT_read_XLSX,mission_actual,conv_UNITS)

Weight_mission = T120_weightconfiguration;
       
switch mission_actual
    %DEPÓSITO LLENO (16.052Kg fuel)
    case 1 % Depósito Lleno 4 misiles
        OUTPUT_read_XLSX.Weights_flags.I_xx = Weight_mission.I_xx_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_yy = Weight_mission.I_yy_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_zz = Weight_mission.I_zz_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_xy = Weight_mission.I_xy_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_xz = Weight_mission.I_xz_vec_f(1); %
        OUTPUT_read_XLSX.Weights_flags.I_yz = Weight_mission.I_yz_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = Weight_mission.x_XCG_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = Weight_mission.y_XCG_vec_f(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = Weight_mission.z_XCG_vec_f(1);
        conditions.m_TOW = Weight_mission.W_vec_f(1);
    case 2 % Depósito Lleno 3 misiles
        OUTPUT_read_XLSX.Weights_flags.I_xx = Weight_mission.I_xx_vec_f(2); %
        OUTPUT_read_XLSX.Weights_flags.I_yy = Weight_mission.I_yy_vec_f(2); %
        OUTPUT_read_XLSX.Weights_flags.I_zz = Weight_mission.I_zz_vec_f(2); %
        OUTPUT_read_XLSX.Weights_flags.I_xy = Weight_mission.I_xy_vec_f(2); %
        OUTPUT_read_XLSX.Weights_flags.I_xz = Weight_mission.I_xz_vec_f(2); %
        OUTPUT_read_XLSX.Weights_flags.I_yz = Weight_mission.I_yz_vec_f(2);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = Weight_mission.x_XCG_vec_f(2);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = Weight_mission.y_XCG_vec_f(2);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = Weight_mission.z_XCG_vec_f(2);
        conditions.m_TOW = Weight_mission.W_vec_f(2);
    case 3 % Depósito Lleno 2 misiles
        OUTPUT_read_XLSX.Weights_flags.I_xx = Weight_mission.I_xx_vec_f(3); %
        OUTPUT_read_XLSX.Weights_flags.I_yy = Weight_mission.I_yy_vec_f(3); %
        OUTPUT_read_XLSX.Weights_flags.I_zz = Weight_mission.I_zz_vec_f(3); %
        OUTPUT_read_XLSX.Weights_flags.I_xy = Weight_mission.I_xy_vec_f(3); %
        OUTPUT_read_XLSX.Weights_flags.I_xz = Weight_mission.I_xz_vec_f(3); %
        OUTPUT_read_XLSX.Weights_flags.I_yz = Weight_mission.I_yz_vec_f(3);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = Weight_mission.x_XCG_vec_f(3);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = Weight_mission.y_XCG_vec_f(3);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = Weight_mission.z_XCG_vec_f(3);
        conditions.m_TOW = Weight_mission.W_vec_f(3);
    case 4 % Depósito Lleno 1 misiles
        OUTPUT_read_XLSX.Weights_flags.I_xx = Weight_mission.I_xx_vec_f(4); %
        OUTPUT_read_XLSX.Weights_flags.I_yy = Weight_mission.I_yy_vec_f(4); %
        OUTPUT_read_XLSX.Weights_flags.I_zz = Weight_mission.I_zz_vec_f(4); %
        OUTPUT_read_XLSX.Weights_flags.I_xy = Weight_mission.I_xy_vec_f(4); %
        OUTPUT_read_XLSX.Weights_flags.I_xz = Weight_mission.I_xz_vec_f(4); %
        OUTPUT_read_XLSX.Weights_flags.I_yz = Weight_mission.I_yz_vec_f(4);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = Weight_mission.x_XCG_vec_f(4);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = Weight_mission.y_XCG_vec_f(4);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = Weight_mission.z_XCG_vec_f(4);
        conditions.m_TOW = Weight_mission.W_vec_f(4);
    case 5 % Depósito Lleno 0 misiles
        OUTPUT_read_XLSX.Weights_flags.I_xx = Weight_mission.I_xx_vec_f(5); %
        OUTPUT_read_XLSX.Weights_flags.I_yy = Weight_mission.I_yy_vec_f(5); %
        OUTPUT_read_XLSX.Weights_flags.I_zz = Weight_mission.I_zz_vec_f(5); %
        OUTPUT_read_XLSX.Weights_flags.I_xy = Weight_mission.I_xy_vec_f(5); %
        OUTPUT_read_XLSX.Weights_flags.I_xz = Weight_mission.I_xz_vec_f(5); %
        OUTPUT_read_XLSX.Weights_flags.I_yz = Weight_mission.I_yz_vec_f(5);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = Weight_mission.x_XCG_vec_f(5);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = Weight_mission.y_XCG_vec_f(5);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = Weight_mission.z_XCG_vec_f(5);
        conditions.m_TOW = Weight_mission.W_vec_f(5);

        % MEDIO DEPÓSITO (8.026KG)
    case 6 % Depósito Lleno 4 misiles
        OUTPUT_read_XLSX.Weights_flags.I_xx = Weight_mission.I_xx_vec_f2(1); %
        OUTPUT_read_XLSX.Weights_flags.I_yy = Weight_mission.I_yy_vec_f2(1); %
        OUTPUT_read_XLSX.Weights_flags.I_zz = Weight_mission.I_zz_vec_f2(1); %
        OUTPUT_read_XLSX.Weights_flags.I_xy = Weight_mission.I_xy_vec_f2(1); %
        OUTPUT_read_XLSX.Weights_flags.I_xz = Weight_mission.I_xz_vec_f2(1); %
        OUTPUT_read_XLSX.Weights_flags.I_yz = Weight_mission.I_yz_vec_f2(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = Weight_mission.x_XCG_vec_f2(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = Weight_mission.y_XCG_vec_f2(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = Weight_mission.z_XCG_vec_f2(1);
        conditions.m_TOW = Weight_mission.W_vec_f2(1);
    case 7 % Depósito Lleno 3 misiles
        OUTPUT_read_XLSX.Weights_flags.I_xx = Weight_mission.I_xx_vec_f2(2); %
        OUTPUT_read_XLSX.Weights_flags.I_yy = Weight_mission.I_yy_vec_f2(2); %
        OUTPUT_read_XLSX.Weights_flags.I_zz = Weight_mission.I_zz_vec_f2(2); %
        OUTPUT_read_XLSX.Weights_flags.I_xy = Weight_mission.I_xy_vec_f2(2); %
        OUTPUT_read_XLSX.Weights_flags.I_xz = Weight_mission.I_xz_vec_f2(2); %
        OUTPUT_read_XLSX.Weights_flags.I_yz = Weight_mission.I_yz_vec_f2(2);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = Weight_mission.x_XCG_vec_f2(2);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = Weight_mission.y_XCG_vec_f2(2);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = Weight_mission.z_XCG_vec_f2(2);
        conditions.m_TOW = Weight_mission.W_vec_f2(2);
    case 8 % Depósito Lleno 2 misiles
        OUTPUT_read_XLSX.Weights_flags.I_xx = Weight_mission.I_xx_vec_f2(3); %
        OUTPUT_read_XLSX.Weights_flags.I_yy = Weight_mission.I_yy_vec_f2(3); %
        OUTPUT_read_XLSX.Weights_flags.I_zz = Weight_mission.I_zz_vec_f2(3); %
        OUTPUT_read_XLSX.Weights_flags.I_xy = Weight_mission.I_xy_vec_f2(3); %
        OUTPUT_read_XLSX.Weights_flags.I_xz = Weight_mission.I_xz_vec_f2(3); %
        OUTPUT_read_XLSX.Weights_flags.I_yz = Weight_mission.I_yz_vec_f2(3);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = Weight_mission.x_XCG_vec_f2(3);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = Weight_mission.y_XCG_vec_f2(3);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = Weight_mission.z_XCG_vec_f2(3);
        conditions.m_TOW = Weight_mission.W_vec_f2(3);
    case 9 % Depósito Lleno 1 misiles
        OUTPUT_read_XLSX.Weights_flags.I_xx = Weight_mission.I_xx_vec_f2(4); %
        OUTPUT_read_XLSX.Weights_flags.I_yy = Weight_mission.I_yy_vec_f2(4); %
        OUTPUT_read_XLSX.Weights_flags.I_zz = Weight_mission.I_zz_vec_f2(4); %
        OUTPUT_read_XLSX.Weights_flags.I_xy = Weight_mission.I_xy_vec_f2(4); %
        OUTPUT_read_XLSX.Weights_flags.I_xz = Weight_mission.I_xz_vec_f2(4); %
        OUTPUT_read_XLSX.Weights_flags.I_yz = Weight_mission.I_yz_vec_f2(4);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = Weight_mission.x_XCG_vec_f2(4);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = Weight_mission.y_XCG_vec_f2(4);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = Weight_mission.z_XCG_vec_f2(4);
        conditions.m_TOW = Weight_mission.W_vec_f2(4);
    case 10 % Depósito Lleno 0 misiles
        OUTPUT_read_XLSX.Weights_flags.I_xx = Weight_mission.I_xx_vec_f2(5); %
        OUTPUT_read_XLSX.Weights_flags.I_yy = Weight_mission.I_yy_vec_f2(5); %
        OUTPUT_read_XLSX.Weights_flags.I_zz = Weight_mission.I_zz_vec_f2(5); %
        OUTPUT_read_XLSX.Weights_flags.I_xy = Weight_mission.I_xy_vec_f2(5); %
        OUTPUT_read_XLSX.Weights_flags.I_xz = Weight_mission.I_xz_vec_f2(5); %
        OUTPUT_read_XLSX.Weights_flags.I_yz = Weight_mission.I_yz_vec_f2(5);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = Weight_mission.x_XCG_vec_f2(5);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = Weight_mission.y_XCG_vec_f2(5);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = Weight_mission.z_XCG_vec_f2(5);
        conditions.m_TOW = Weight_mission.W_vec_f2(5);

        % DEPÓSITO VACÍO
    case 11 % Depósito Lleno 4 misiles
        OUTPUT_read_XLSX.Weights_flags.I_xx = Weight_mission.I_xx_vec_nf(1); %
        OUTPUT_read_XLSX.Weights_flags.I_yy = Weight_mission.I_yy_vec_nf(1); %
        OUTPUT_read_XLSX.Weights_flags.I_zz = Weight_mission.I_zz_vec_nf(1); %
        OUTPUT_read_XLSX.Weights_flags.I_xy = Weight_mission.I_xy_vec_nf(1); %
        OUTPUT_read_XLSX.Weights_flags.I_xz = Weight_mission.I_xz_vec_nf(1); %
        OUTPUT_read_XLSX.Weights_flags.I_yz = Weight_mission.I_yz_vec_nf(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = Weight_mission.x_XCG_vec_nf(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = Weight_mission.y_XCG_vec_nf(1);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = Weight_mission.z_XCG_vec_nf(1);
        conditions.m_TOW = Weight_mission.W_vec_nf(1);
    case 12 % Depósito Lleno 3 misiles
        OUTPUT_read_XLSX.Weights_flags.I_xx = Weight_mission.I_xx_vec_nf(2); %
        OUTPUT_read_XLSX.Weights_flags.I_yy = Weight_mission.I_yy_vec_nf(2); %
        OUTPUT_read_XLSX.Weights_flags.I_zz = Weight_mission.I_zz_vec_nf(2); %
        OUTPUT_read_XLSX.Weights_flags.I_xy = Weight_mission.I_xy_vec_nf(2); %
        OUTPUT_read_XLSX.Weights_flags.I_xz = Weight_mission.I_xz_vec_nf(2); %
        OUTPUT_read_XLSX.Weights_flags.I_yz = Weight_mission.I_yz_vec_nf(2);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = Weight_mission.x_XCG_vec_nf(2);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = Weight_mission.y_XCG_vec_nf(2);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = Weight_mission.z_XCG_vec_nf(2);
        conditions.m_TOW = Weight_mission.W_vec_nf(2);
    case 13 % Depósito Lleno 2 misiles
        OUTPUT_read_XLSX.Weights_flags.I_xx = Weight_mission.I_xx_vec_nf(3); %
        OUTPUT_read_XLSX.Weights_flags.I_yy = Weight_mission.I_yy_vec_nf(3); %
        OUTPUT_read_XLSX.Weights_flags.I_zz = Weight_mission.I_zz_vec_nf(3); %
        OUTPUT_read_XLSX.Weights_flags.I_xy = Weight_mission.I_xy_vec_nf(3); %
        OUTPUT_read_XLSX.Weights_flags.I_xz = Weight_mission.I_xz_vec_nf(3); %
        OUTPUT_read_XLSX.Weights_flags.I_yz = Weight_mission.I_yz_vec_nf(3);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = Weight_mission.x_XCG_vec_nf(3);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = Weight_mission.y_XCG_vec_nf(3);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = Weight_mission.z_XCG_vec_nf(3);
        conditions.m_TOW = Weight_mission.W_vec_nf(3);
    case 14 % Depósito Lleno 1 misiles
        OUTPUT_read_XLSX.Weights_flags.I_xx = Weight_mission.I_xx_vec_nf(4); %
        OUTPUT_read_XLSX.Weights_flags.I_yy = Weight_mission.I_yy_vec_nf(4); %
        OUTPUT_read_XLSX.Weights_flags.I_zz = Weight_mission.I_zz_vec_nf(4); %
        OUTPUT_read_XLSX.Weights_flags.I_xy = Weight_mission.I_xy_vec_nf(4); %
        OUTPUT_read_XLSX.Weights_flags.I_xz = Weight_mission.I_xz_vec_nf(4); %
        OUTPUT_read_XLSX.Weights_flags.I_yz = Weight_mission.I_yz_vec_nf(4);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = Weight_mission.x_XCG_vec_nf(4);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = Weight_mission.y_XCG_vec_nf(4);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = Weight_mission.z_XCG_vec_nf(4);
        conditions.m_TOW = Weight_mission.W_vec_nf(4);
    case 15 % Depósito Lleno 0 misiles
        OUTPUT_read_XLSX.Weights_flags.I_xx = Weight_mission.I_xx_vec_nf(5); %
        OUTPUT_read_XLSX.Weights_flags.I_yy = Weight_mission.I_yy_vec_nf(5); %
        OUTPUT_read_XLSX.Weights_flags.I_zz = Weight_mission.I_zz_vec_nf(5); %
        OUTPUT_read_XLSX.Weights_flags.I_xy = Weight_mission.I_xy_vec_nf(5); %
        OUTPUT_read_XLSX.Weights_flags.I_xz = Weight_mission.I_xz_vec_nf(5); %
        OUTPUT_read_XLSX.Weights_flags.I_yz = Weight_mission.I_yz_vec_nf(5);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.x_XCG = Weight_mission.x_XCG_vec_nf(5);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.y_XCG = Weight_mission.y_XCG_vec_nf(5);
        OUTPUT_read_XLSX.InputGeometry_Data_flags.z_XCG = Weight_mission.z_XCG_vec_nf(5);
        conditions.m_TOW = Weight_mission.W_vec_nf(5);
end
