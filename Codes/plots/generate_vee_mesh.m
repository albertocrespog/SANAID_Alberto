function [PLOTTING_UAV,x_vee_LE] = generate_vee_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX,AC_CONFIGURATION)

W1 = AC_CONFIGURATION.W1;
HTP = AC_CONFIGURATION.HTP;
VTP = AC_CONFIGURATION.VTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
twin_VTP = AC_CONFIGURATION.twin_VTP;

% Surfaces
S_vee = Geo_tier.S_vee;

% Spans
b_vee = Geo_tier.b_vee;
b_vee_e = Geo_tier.b_vee_e;

% Root Chords
cR_vee = Geo_tier.cR_vee;

% Tip Chords
cT_vee = Geo_tier.cT_vee;

% mean chord
c_vee = Geo_tier.cmac_vee;

dihedral_vee = Geo_tier.dihedral_vee;

% Sweep
Lambda_LE_vee = Geo_tier.Lambda_LE_vee;
Lambda_TE_vee = Geo_tier.Lambda_TE_vee;
Lambda_c4_vee = Geo_tier.Lambda_c4_vee;
Lambda_c2_vee = Geo_tier.Lambda_c2_vee;

% %------------------------------GEOMETRY------------------------------%
% % posición del borde de ataque de las superficies aerodinámicas
x_vee_LE = Geo_tier.x_vee_LE;

y_offset_vee = Geo_tier.y_offset_vee;

x_cR_vee_LE = Geo_tier.x_cR_vee_LE;
x_cR_vee_TE = Geo_tier.x_cR_vee_TE;
x_cT_vee_LE = Geo_tier.x_cT_vee_LE;
x_cT_vee_TE = Geo_tier.x_cT_vee_TE;

y_cR_vee_LE = Geo_tier.y_cR_vee_LE;
z_cR_vee_LE = Geo_tier.z_cR_vee_LE;
z_cR_vee_TE = Geo_tier.z_cR_vee_TE;
z_cT_vee_LE = Geo_tier.z_cT_vee_LE;

%--------------------------- WING 2 ---------------------------------
S = S_vee;
b = b_vee/2;
b_vee_fus_in =(0); % increase in y-direction of span associated to fuselage width (root)
b_vee_fus_out = (b_vee_e); % increase in y-direction of span associated to fuselage width (tip)
y = [b_vee_fus_in b_vee_fus_out]'/2;
c_y = [cR_vee cT_vee]';
le_y_1 = (x_cR_vee_LE - x_cR_vee_TE);
le_y_2 = (x_cT_vee_LE - x_cR_vee_TE);
le_y = [le_y_1 le_y_2]';
diedro = [dihedral_vee dihedral_vee];
NACA_foil = 1;
VTP_ms = 0;
t_c = 0.10  ;
[x_mesh_vee, y_mesh_vee, z_mesh_vee,  CA_vee, MISC_vee] = get_DATA_aero(NACA_foil,S,b,y,c_y,le_y,diedro,VTP_ms,t_c);

% Wing2
X = x_cR_vee_TE;
Y = y_offset_vee;
center_section = 1;
Z = z_cR_vee_TE;
[x_mesh_vee_New, y_mesh_vee_New, z_mesh_vee_New] = get_DATA_New_Location(x_mesh_vee, y_mesh_vee, z_mesh_vee,X,Y,Z,center_section);


PLOTTING_UAV.x_mesh_vee_New = x_mesh_vee_New; % vee
PLOTTING_UAV.y_mesh_vee_New = y_mesh_vee_New; % vee
PLOTTING_UAV.z_mesh_vee_New = z_mesh_vee_New; % vee
