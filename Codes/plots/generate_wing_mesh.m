function [PLOTTING_UAV,x_w1_LE] = generate_wing_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX,AC_CONFIGURATION)

W1 = AC_CONFIGURATION.W1;
HTP = AC_CONFIGURATION.HTP;
VTP = AC_CONFIGURATION.VTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
twin_VTP = AC_CONFIGURATION.twin_VTP;

% Surfaces
S_w1 = Geo_tier.S_w1;

% Spans
b_w1 = Geo_tier.b_w1;
b_w1_e = Geo_tier.b_w1_e;

% Root Chords
cR_w1 = Geo_tier.cR_w1;
% Tip Chords
cT_w1 = Geo_tier.cT_w1;

% mean chord
c_w1 = Geo_tier.cmac_w1;

dihedral_w1 = Geo_tier.dihedral_w1;

% Sweep
Lambda_LE_w1 = Geo_tier.Lambda_LE_w1;
Lambda_TE_w1 = Geo_tier.Lambda_TE_w1;
Lambda_c4_w1 = Geo_tier.Lambda_c4_w1;
Lambda_c2_w1 = Geo_tier.Lambda_c2_w1;

% if OUTPUT_read_XLSX.Fuselage_flags.CAD_Kink_Wing == 1
%     % Wingspan
%     y_loc_1R_y1_w1_CAD = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_y1_w1_CAD;
%     y_loc_1R_yB1_w1_CAD = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_yB1_w1_CAD;
%     y_loc_1R_yB2_w1_CAD = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_loc_1R_yB2_w1_CAD;
%     % Sweep
%     Lambda_LE_w1_e = OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w1_e;
%     Lambda_LE_w1_k1_e = OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w1_k1_e;
%     Lambda_LE_w1_k2_e = OUTPUT_read_XLSX.InputGeometry_Data_flags.Lambda_LE_w1_k2_e;
%     % Dihedral
%     dihedral_w1_e = OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w1_e;
%     dihedral_w1_k1_e = OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w1_k1_e;
%     dihedral_w1_k2_e = OUTPUT_read_XLSX.InputGeometry_Data_flags.dihedral_w1_k2_e;
%     % Chrod
%     cR_w1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.cR_w1;
%     cB_k1_w1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.cB_k1_w1;
%     cB_k2_w1 = OUTPUT_read_XLSX.InputGeometry_Data_flags.cB_k2_w1;
% 
%     % %------------------------------GEOMETRY------------------------------%
%     % % posición del borde de ataque de las superficies aerodinámicas
%     x_w1_LE = Geo_tier.x_w1_LE;
% 
%     y_offset_w1 = Geo_tier.y_offset_w1;
% 
%     % KINK1
%     x_cR_w1_LE = Geo_tier.x_cR_w1_LE;
%     x_cR_w1_TE = Geo_tier.x_cR_w1_TE;
%     y_cR_w1_LE = Geo_tier.y_cR_w1_LE;
%     y_cR_w1_TE = Geo_tier.y_cR_w1_TE;
%     z_cR_w1_LE = Geo_tier.z_cR_w1_LE;
%     z_cR_w1_TE = Geo_tier.z_cR_w1_TE;
% 
%     x_cT_w1_LE = Geo_tier.x_cT_w1_LE;
%     x_cT_w1_TE = Geo_tier.x_cT_w1_TE;
%     y_cT_w1_LE = Geo_tier.y_cT_w1_LE;
%     y_cT_w1_TE = Geo_tier.y_cT_w1_TE;
%     z_cT_w1_LE = Geo_tier.z_cT_w1_LE;
%     z_cT_w1_TE = Geo_tier.z_cT_w1_TE;
% 
% else

    % %------------------------------GEOMETRY------------------------------%
    % % posición del borde de ataque de las superficies aerodinámicas
    x_w1_LE = Geo_tier.x_w1_LE;

    y_offset_w1 = Geo_tier.y_offset_w1;

    x_cR_w1_LE = Geo_tier.x_cR_w1_LE;
    x_cR_w1_TE = Geo_tier.x_cR_w1_TE;
    y_cR_w1_LE = Geo_tier.y_cR_w1_LE;
    y_cR_w1_TE = Geo_tier.y_cR_w1_TE;
    z_cR_w1_LE = Geo_tier.z_cR_w1_LE;
    z_cR_w1_TE = Geo_tier.z_cR_w1_TE;

    x_cT_w1_LE = Geo_tier.x_cT_w1_LE;
    x_cT_w1_TE = Geo_tier.x_cT_w1_TE;
    y_cT_w1_LE = Geo_tier.y_cT_w1_LE;
    y_cT_w1_TE = Geo_tier.y_cT_w1_TE;
    z_cT_w1_LE = Geo_tier.z_cT_w1_LE;
    z_cT_w1_TE = Geo_tier.z_cT_w1_TE;

    %--------------------------- WING 1 ---------------------------------
    S = S_w1;
    b = b_w1/2;
    b_w1_fus_in =(0); % increase in y-direction of span associated to fuselage width (root)
    b_w1_fus_out = (b_w1_e); % increase in y-direction of span associated to fuselage width (tip)
    y = [b_w1_fus_in b_w1_fus_out]'/2;
    c_y = [cR_w1 cT_w1]';
    le_y_1 = (x_cR_w1_LE - x_cR_w1_TE);
    le_y_2 = (x_cT_w1_LE - x_cR_w1_TE);
    le_y = [le_y_1 le_y_2]';
    diedro = [dihedral_w1 dihedral_w1];
    NACA_foil = 1;
    VTP_ms = 0;
    t_c = 0.18;
    [x_mesh_w1,y_mesh_w1,z_mesh_w1,CA_w1,MISC_w1] = get_DATA_aero(NACA_foil,S,b,y,c_y,le_y,diedro,VTP_ms,t_c);

    % Modifies location of Mesh
    % Wing1
    X = x_cR_w1_TE;
    Y = y_offset_w1;
    Y = 0;
    Y=y_cR_w1_LE;
    center_section = 1; % Flag that determine if symmetry along the x axis needs to be appliesd
    Z = z_cR_w1_TE;
    [x_mesh_w1_New, y_mesh_w1_New, z_mesh_w1_New] = get_DATA_New_Location(x_mesh_w1, y_mesh_w1, z_mesh_w1,X,Y,Z,center_section);

    PLOTTING_UAV.x_mesh_w1_New = x_mesh_w1_New; % WING
    PLOTTING_UAV.y_mesh_w1_New = y_mesh_w1_New; % WING
    PLOTTING_UAV.z_mesh_w1_New = z_mesh_w1_New; % WING
% end
