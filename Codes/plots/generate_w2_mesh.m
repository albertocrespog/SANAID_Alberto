function [PLOTTING_UAV,x_w2_LE] = generate_w2_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX,AC_CONFIGURATION)

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
Lambda_LE_w2 = Geo_tier.Lambda_LE_w2;
Lambda_TE_w2 = Geo_tier.Lambda_TE_w2;
Lambda_c4_w2 = Geo_tier.Lambda_c4_w2;
Lambda_c2_w2 = Geo_tier.Lambda_c2_w2;

% %------------------------------GEOMETRY------------------------------%
% % posición del borde de ataque de las superficies aerodinámicas
x_w2_LE = Geo_tier.x_w2_LE;

y_offset_w2 = Geo_tier.y_offset_w2;

x_cR_w2_LE = Geo_tier.x_cR_w2_LE;
x_cR_w2_TE = Geo_tier.x_cR_w2_TE;
x_cT_w2_LE = Geo_tier.x_cT_w2_LE;
x_cT_w2_TE = Geo_tier.x_cT_w2_TE;

y_cR_w2_LE = Geo_tier.y_cR_w2_LE;
z_cR_w2_LE = Geo_tier.z_cR_w2_LE;
z_cR_w2_TE = Geo_tier.z_cR_w2_TE;
z_cT_w2_LE = Geo_tier.z_cT_w2_LE;

%--------------------------- WING 2 ---------------------------------
S = S_w2;
b = b_w2/2;
b_w2_fus_in =(0); % increase in y-direction of span associated to fuselage width (root)
b_w2_fus_out = (b_w2_e); % increase in y-direction of span associated to fuselage width (tip)
y = [b_w2_fus_in b_w2_fus_out]'/2;
c_y = [cR_w2 cT_w2]';
le_y_1 = (x_cR_w2_LE - x_cR_w2_TE);
le_y_2 = (x_cT_w2_LE - x_cR_w2_TE);
le_y = [le_y_1 le_y_2]';
diedro = [dihedral_w2 dihedral_w2];
NACA_foil = 1;
VTP_ms = 0;
t_c = 0.10  ;
[x_mesh_w2, y_mesh_w2, z_mesh_w2,  CA_w2, MISC_w2] = get_DATA_aero(NACA_foil,S,b,y,c_y,le_y,diedro,VTP_ms,t_c);

% Wing2
X = x_cR_w2_TE;
Y = y_offset_w2;
center_section = 1;
Z = z_cR_w2_TE;
[x_mesh_w2_New, y_mesh_w2_New, z_mesh_w2_New] = get_DATA_New_Location(x_mesh_w2, y_mesh_w2, z_mesh_w2,X,Y,Z,center_section);


PLOTTING_UAV.x_mesh_w2_New = x_mesh_w2_New; % HTP
PLOTTING_UAV.y_mesh_w2_New = y_mesh_w2_New; % HTP
PLOTTING_UAV.z_mesh_w2_New = z_mesh_w2_New; % HTP
