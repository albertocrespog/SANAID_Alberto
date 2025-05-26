function [PLOTTING_UAV,x_vee2_LE] = generate_vee2_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX,AC_CONFIGURATION)

W1 = AC_CONFIGURATION.W1;
HTP = AC_CONFIGURATION.HTP;
VTP = AC_CONFIGURATION.VTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
twin_VTP = AC_CONFIGURATION.twin_VTP;

% Surfaces
S_vee2 = Geo_tier.S_vee2;

% Spans
b_vee2 = Geo_tier.b_vee;
b_vee2_e = Geo_tier.b_vee2_e;

% Root Chords
cR_vee2 = Geo_tier.cR_vee2;

% Tip Chords
cT_vee2 = Geo_tier.cT_vee2;

% mean chord
c_vee2 = Geo_tier.cmac_vee2;

dihedral_vee2 = Geo_tier.dihedral_vee2;

% Sweep
Lambda_LE_vee2 = Geo_tier.Lambda_LE_vee2;
Lambda_TE_vee2 = Geo_tier.Lambda_TE_vee2;
Lambda_c4_vee2 = Geo_tier.Lambda_c4_vee2;
Lambda_c2_vee2 = Geo_tier.Lambda_c2_vee2;

% %------------------------------GEOMETRY------------------------------%
% % posición del borde de ataque de las superficies aerodinámicas
x_vee2_LE = Geo_tier.x_vee2_LE;

y_offset_vee2 = Geo_tier.y_offset_vee2;

x_cR_vee2_LE = Geo_tier.x_cR_vee2_LE;
x_cR_vee2_TE = Geo_tier.x_cR_vee2_TE;
x_cT_vee2_LE = Geo_tier.x_cT_vee2_LE;
x_cT_vee2_TE = Geo_tier.x_cT_vee2_TE;

y_cR_vee2_LE = Geo_tier.y_cR_vee2_LE;
z_cR_vee2_LE = Geo_tier.z_cR_vee2_LE;
z_cR_vee2_TE = Geo_tier.z_cR_vee2_TE;
z_cT_vee2_LE = Geo_tier.z_cT_vee2_LE;

%--------------------------- WING 2 ---------------------------------
S = S_vee2;
b = b_vee2/2;
b_vee2_fus_in =(0); % increase in y-direction of span associated to fuselage width (root)
b_vee2_fus_out = (b_vee2_e); % increase in y-direction of span associated to fuselage width (tip)
y = [b_vee2_fus_in b_vee2_fus_out]'/2;
c_y = [cR_vee2 cT_vee2]';
le_y_1 = (x_cR_vee2_LE - x_cR_vee2_TE);
le_y_2 = (x_cT_vee2_LE - x_cR_vee2_TE);
le_y = [le_y_1 le_y_2]';
diedro = [dihedral_vee2 dihedral_vee2];
NACA_foil = 1;
VTP_ms = 0;
t_c = 0.10  ;
[x_mesh_vee2, y_mesh_vee2, z_mesh_vee2,  CA_vee2, MISC_vee2] = get_DATA_aero(NACA_foil,S,b,y,c_y,le_y,diedro,VTP_ms,t_c);

% Wing2
X = x_cR_vee2_TE;
Y = y_offset_vee2;
center_section = 1;
Z = z_cR_vee2_TE;
[x_mesh_vee2_New, y_mesh_vee2_New, z_mesh_vee2_New] = get_DATA_New_Location(x_mesh_vee2, y_mesh_vee2, z_mesh_vee2,X,Y,Z,center_section);


PLOTTING_UAV.x_mesh_vee2_New = x_mesh_vee2_New; % vee2
PLOTTING_UAV.y_mesh_vee2_New = y_mesh_vee2_New; % vee2
PLOTTING_UAV.z_mesh_vee2_New = z_mesh_vee2_New; % vee2
