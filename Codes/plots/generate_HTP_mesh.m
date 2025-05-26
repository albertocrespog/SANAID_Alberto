function [PLOTTING_UAV,x_HTP_LE] = generate_HTP_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX,AC_CONFIGURATION)

W1 = AC_CONFIGURATION.W1;
HTP = AC_CONFIGURATION.HTP;
VTP = AC_CONFIGURATION.VTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
twin_VTP = AC_CONFIGURATION.twin_VTP;

% Surfaces
S_HTP = Geo_tier.S_HTP;

% Spans
b_HTP = Geo_tier.b_HTP;
b_HTP_e = Geo_tier.b_HTP_e;

% Root Chords
cR_HTP = Geo_tier.cR_HTP;

% Tip Chords
cT_HTP = Geo_tier.cT_HTP;

% mean chord
c_HTP = Geo_tier.cmac_HTP;

dihedral_HTP = Geo_tier.dihedral_HTP;

% Sweep
Lambda_LE_HTP = Geo_tier.Lambda_LE_HTP;
Lambda_TE_HTP = Geo_tier.Lambda_TE_HTP;
Lambda_c4_HTP = Geo_tier.Lambda_c4_HTP;
Lambda_c2_HTP = Geo_tier.Lambda_c2_HTP;

% %------------------------------GEOMETRY------------------------------%
% % posición del borde de ataque de las superficies aerodinámicas
x_HTP_LE = Geo_tier.x_HTP_LE;

y_offset_HTP = Geo_tier.y_offset_HTP;

x_cR_HTP_LE = Geo_tier.x_cR_HTP_LE;
x_cR_HTP_TE = Geo_tier.x_cR_HTP_TE;
x_cT_HTP_LE = Geo_tier.x_cT_HTP_LE;
x_cT_HTP_TE = Geo_tier.x_cT_HTP_TE;

y_cR_HTP_LE = Geo_tier.y_cR_HTP_LE;
z_cR_HTP_LE = Geo_tier.z_cR_HTP_LE;
z_cR_HTP_TE = Geo_tier.z_cR_HTP_TE;
z_cT_HTP_LE = Geo_tier.z_cT_HTP_LE;

%--------------------------- WING 2 ---------------------------------
S = S_HTP;
b = b_HTP/2;
b_HTP_fus_in =(0); % increase in y-direction of span associated to fuselage width (root)
b_HTP_fus_out = (b_HTP_e); % increase in y-direction of span associated to fuselage width (tip)
y = [b_HTP_fus_in b_HTP_fus_out]'/2;
c_y = [cR_HTP cT_HTP]';
le_y_1 = (x_cR_HTP_LE - x_cR_HTP_TE);
le_y_2 = (x_cT_HTP_LE - x_cR_HTP_TE);
le_y = [le_y_1 le_y_2]';
diedro = [dihedral_HTP dihedral_HTP];
NACA_foil = 1;
VTP_ms = 0;
t_c = 0.10  ;
[x_mesh_HTP, y_mesh_HTP, z_mesh_HTP,  CA_HTP, MISC_HTP] = get_DATA_aero(NACA_foil,S,b,y,c_y,le_y,diedro,VTP_ms,t_c);

% Wing2
X = x_cR_HTP_TE;
Y = y_offset_HTP;
center_section = 1;
Z = z_cR_HTP_TE;
[x_mesh_HTP_New, y_mesh_HTP_New, z_mesh_HTP_New] = get_DATA_New_Location(x_mesh_HTP, y_mesh_HTP, z_mesh_HTP,X,Y,Z,center_section);


PLOTTING_UAV.x_mesh_HTP_New = x_mesh_HTP_New; % HTP
PLOTTING_UAV.y_mesh_HTP_New = y_mesh_HTP_New; % HTP
PLOTTING_UAV.z_mesh_HTP_New = z_mesh_HTP_New; % HTP
