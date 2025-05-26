function [PLOTTING_UAV, x_can_LE] = generate_can_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX,AC_CONFIGURATION)

W1 = AC_CONFIGURATION.W1;
HTP = AC_CONFIGURATION.HTP;
VTP = AC_CONFIGURATION.VTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
twin_VTP = AC_CONFIGURATION.twin_VTP;

% Surfaces
S_can = Geo_tier.S_can;

% Spans
b_can = Geo_tier.b_can;
b_can_e = Geo_tier.b_can_e;

% Root Chords
cR_can = Geo_tier.cR_can;

% Tip Chords
cT_can = Geo_tier.cT_can;

% mean chord
c_can = Geo_tier.cmac_can;

dihedral_can = Geo_tier.dihedral_can;

% Sweep
Lambda_LE_can = Geo_tier.Lambda_LE_can;
Lambda_TE_can = Geo_tier.Lambda_TE_can;
Lambda_c4_can = Geo_tier.Lambda_c4_can;
Lambda_c2_can = Geo_tier.Lambda_c2_can;

% %------------------------------GEOMETRY------------------------------%
% % posición del borde de ataque de las superficies aerodinámicas
x_can_LE = Geo_tier.x_can_LE;

y_offset_can = Geo_tier.y_offset_can;

x_cR_can_LE = Geo_tier.x_cR_can_LE;
x_cR_can_TE = Geo_tier.x_cR_can_TE;
x_cT_can_LE = Geo_tier.x_cT_can_LE;
x_cT_can_TE = Geo_tier.x_cT_can_TE;

y_cR_can_LE = Geo_tier.y_cR_can_LE;
z_cR_can_LE = Geo_tier.z_cR_can_LE;
z_cR_can_TE = Geo_tier.z_cR_can_TE;
z_cT_can_LE = Geo_tier.z_cT_can_LE;

%--------------------------- CANARD ---------------------------------
S = S_can;
b = b_can/2;
b_can_fus_in =(0); % increase in y-direction of span associated to fuselage width (root)
b_can_fus_out = (b_can_e); % increase in y-direction of span associated to fuselage width (tip)
y = [b_can_fus_in b_can_fus_out]'/2;
c_y = [cR_can cT_can]';
le_y_1 = (x_cR_can_LE - x_cR_can_TE);
le_y_2 = (x_cT_can_LE - x_cR_can_TE);
le_y = [le_y_1 le_y_2]';
diedro = [dihedral_can dihedral_can];
NACA_foil = 1;
VTP_ms = 0;
t_c = 0.18;
[x_mesh_can,y_mesh_can,z_mesh_can,CA_can,MISC_can] = get_DATA_aero(NACA_foil,S,b,y,c_y,le_y,diedro,VTP_ms,t_c);

% Modifies location of Mesh
% Wing1
X = x_cR_can_TE;
if OUTPUT_read_XLSX.InputGeometry_Data_flags.canard_offset_can == 1
    Y = y_offset_can;
else
    Y = 0;
end
center_section = 1; % Flag that determine if symmetry along the x axis needs to be appliesd
Z = z_cR_can_TE;
[x_mesh_can_New, y_mesh_can_New, z_mesh_can_New] = get_DATA_New_Location(x_mesh_can, y_mesh_can, z_mesh_can,X,Y,Z,center_section);

PLOTTING_UAV.x_mesh_can_New = x_mesh_can_New; % Canard
PLOTTING_UAV.y_mesh_can_New = y_mesh_can_New; % Canard
PLOTTING_UAV.z_mesh_can_New = z_mesh_can_New; % Canard
