function [PLOTTING_UAV, x_VTP_LE] = generate_VTP_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX,AC_CONFIGURATION)

W1 = AC_CONFIGURATION.W1;
HTP = AC_CONFIGURATION.HTP;
VTP = AC_CONFIGURATION.VTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
twin_VTP = AC_CONFIGURATION.twin_VTP;

% Surfaces
S_VTP = Geo_tier.S_VTP;

% Spans
b_VTP = Geo_tier.b_VTP;
b_VTP_e = Geo_tier.b_VTP_e;

% Root Chords
cR_VTP = Geo_tier.cR_VTP;

% Tip Chords
cT_VTP = Geo_tier.cT_VTP;

% mean chord
c_VTP = Geo_tier.cmac_VTP;

dihedral_VTP = Geo_tier.dihedral_VTP;

% Sweep
Lambda_LE_VTP = Geo_tier.Lambda_LE_VTP;
Lambda_TE_VTP = Geo_tier.Lambda_TE_VTP;
Lambda_c4_VTP = Geo_tier.Lambda_c4_VTP;
Lambda_c2_VTP = Geo_tier.Lambda_c2_VTP;

% %------------------------------GEOMETRY------------------------------%
% % posición del borde de ataque de las superficies aerodinámicas
x_VTP_LE = Geo_tier.x_VTP_LE;

y_offset_VTP = Geo_tier.y_offset_VTP;

x_cR_VTP_LE = Geo_tier.x_cR_VTP_LE;
x_cR_VTP_TE = Geo_tier.x_cR_VTP_TE;
y_cR_VTP_LE = Geo_tier.y_cR_VTP_LE;
y_cR_VTP_TE = Geo_tier.y_cR_VTP_TE;
z_cR_VTP_LE = Geo_tier.z_cR_VTP_LE;
z_cR_VTP_TE = Geo_tier.z_cR_VTP_TE;

x_cT_VTP_LE = Geo_tier.x_cT_VTP_LE;
x_cT_VTP_TE = Geo_tier.x_cT_VTP_TE;
y_cT_VTP_LE = Geo_tier.y_cT_VTP_LE;
y_cT_VTP_TE = Geo_tier.y_cT_VTP_TE;
z_cT_VTP_LE = Geo_tier.z_cT_VTP_LE;
z_cT_VTP_TE = Geo_tier.z_cT_VTP_TE;

%--------------------------- VTP ---------------------------------
S = S_VTP;
b = b_VTP;
b_VTP_fus_in =(0); % increase in y-direction of span associated to fuselage width (root)
b_VTP_fus_out = (b_VTP_e); % increase in y-direction of span associated to fuselage width (tip)
y = [b_VTP_fus_in b_VTP_fus_out]';
c_y = [cR_VTP cT_VTP]';
le_y_1 = (x_cR_VTP_LE - x_cR_VTP_TE);
le_y_2 = (x_cT_VTP_LE - x_cR_VTP_TE);
le_y = [le_y_1 le_y_2]';
dihedral_VTP = 0;
diedro = [dihedral_VTP dihedral_VTP];
NACA_foil = 1;
VTP_ms = 1;
t_c = 0.10;
[x_mesh_VTP_tmp, y_mesh_VTP_tmp, z_mesh_VTP_tmp,  CA_VTP, MISC_VTP] = get_DATA_aero(NACA_foil,S,b,y,c_y,le_y,diedro,VTP_ms,t_c);

x_mesh_VTP =  x_mesh_VTP_tmp;
y_mesh_VTP =  y_mesh_VTP_tmp;
z_mesh_VTP =  z_mesh_VTP_tmp;

% Modifies location of Mesh
if VTP == 1
    if twin_VTP == 1

        % VTP1
        center_section = 0;
        X = x_cR_VTP_TE;
        Y = y_offset_VTP;
        Z = z_cR_VTP_TE;
        [x_mesh_VTP_New_tmp, y_mesh_VTP_New_tmp, z_mesh_VTP_New_tmp] = get_DATA_New_Location(x_mesh_VTP, z_mesh_VTP, y_mesh_VTP,X,Y,Z,center_section);

        x_mesh_VTP1_New =  x_mesh_VTP_New_tmp;
        y_mesh_VTP1_New =  y_mesh_VTP_New_tmp;
        z_mesh_VTP1_New =  z_mesh_VTP_New_tmp;

        % VTP2
        center_section = 0;
        X = x_cR_VTP_TE;
        Y = -y_offset_VTP;
        Z = z_cR_VTP_TE;
        [x_mesh_VTP_New_tmp, y_mesh_VTP_New_tmp, z_mesh_VTP_New_tmp] = get_DATA_New_Location(x_mesh_VTP, z_mesh_VTP, y_mesh_VTP,X,Y,Z,center_section);

        x_mesh_VTP2_New =  x_mesh_VTP_New_tmp;
        y_mesh_VTP2_New =  y_mesh_VTP_New_tmp;
        z_mesh_VTP2_New =  z_mesh_VTP_New_tmp;
        
        PLOTTING_UAV.x_mesh_VTP1_New = x_mesh_VTP1_New; % VTP1
        PLOTTING_UAV.y_mesh_VTP1_New = y_mesh_VTP1_New; % VTP1
        PLOTTING_UAV.z_mesh_VTP1_New = z_mesh_VTP1_New; % VTP1
        PLOTTING_UAV.x_mesh_VTP2_New = x_mesh_VTP2_New; % VTP2
        PLOTTING_UAV.y_mesh_VTP2_New = y_mesh_VTP2_New; % VTP2
        PLOTTING_UAV.z_mesh_VTP2_New = z_mesh_VTP2_New; % VTP2

    else
        center_section = 0;
        X = x_cR_VTP_TE;
        Y = y_offset_VTP;
        Z = z_cR_VTP_TE;
        [x_mesh_VTP_New_tmp, y_mesh_VTP_New_tmp, z_mesh_VTP_New_tmp] = get_DATA_New_Location(x_mesh_VTP,z_mesh_VTP,y_mesh_VTP,X,Y,Z,center_section);

        x_mesh_VTP_New =  x_mesh_VTP_New_tmp;
        y_mesh_VTP_New =  y_mesh_VTP_New_tmp;
        z_mesh_VTP_New =  z_mesh_VTP_New_tmp;

        PLOTTING_UAV.x_mesh_VTP_New = x_mesh_VTP_New; % VTP
        PLOTTING_UAV.y_mesh_VTP_New = y_mesh_VTP_New; % VTP
        PLOTTING_UAV.z_mesh_VTP_New = z_mesh_VTP_New; % VTP
    end
end