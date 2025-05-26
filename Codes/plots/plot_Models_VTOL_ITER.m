function PLOTTING_UAV = plot_Models_VTOL_ITER(Geo_tier,Body_Geo,meshData,Prop_data,AC_CONFIGURATION,OUTPUT_read_XLSX)

% Constants
f2m = 0.3048;
D2R = pi/180;
R2D = 180/pi;

n_eng = Prop_data.n_eng;

% Stores Aircraft Configuration
% AC_CONFIGURATION.Control_surface = Control_surface;
AC_type = AC_CONFIGURATION.AC_type;
Engine_loc = AC_CONFIGURATION.Engine_loc;
Engine_conf = AC_CONFIGURATION.Engine_conf;

W1 = AC_CONFIGURATION.W1;
HTP = AC_CONFIGURATION.HTP;
VTP = AC_CONFIGURATION.VTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
twin_VTP = AC_CONFIGURATION.twin_VTP;

% Fuselage geometry
D_fus = Geo_tier.d_fus;
W_fus = Geo_tier.w_fus;
H_fus = Geo_tier.h_fus;
L_fus = Geo_tier.l_fus;

% Dummy so that it overwrites the variables
PLOTTING_UAV.dummy = 0;

% Aircraft type
switch AC_type
    case 1 % AC_type = 1 - flying wing

%         generate_wing_mesh(Geo_tier)
% 
%         % Surfaces
%         S_w1 = Geo_tier.S_w1;
%         
%         % Spans
%         b_w1 = Geo_tier.b_w1;
%         b_w1_e = Geo_tier.b_w1_e;
%         
%         % Root Chords
%         cR_w1 = Geo_tier.cR_w1;
%         cR_w2 = Geo_tier.cR_w2;
%         
%         % Tip Chords
%         cT_w1 = Geo_tier.cT_w1;
%         
%         % mean chord
%         c_w1 = Geo_tier.cmac_w1;
%         c_w2 = Geo_tier.cmac_w2;
%         
%         dihedral_w1 = Geo_tier.dihedral_w1;
% 
%         % Sweep
%         Lambda_LE_w1 = Geo_tier.Lambda_LE_w1;
%         Lambda_TE_w1 = Geo_tier.Lambda_TE_w1;
%         Lambda_c4_w1 = Geo_tier.Lambda_c4_w1;
%         Lambda_c2_w1 = Geo_tier.Lambda_c2_w1;
%         
%         % %------------------------------GEOMETRY------------------------------%
%         % % posición del borde de ataque de las superficies aerodinámicas
%         x_w1_LE = Geo_tier.x_w1_LE;
%         
%         y_offset_w1 = Geo_tier.y_offset_w1;
%         
%         x_cR_w1_LE = Geo_tier.x_cR_w1_LE;
%         x_cR_w1_TE = Geo_tier.x_cR_w1_TE;
%         x_cT_w1_LE = Geo_tier.x_cT_w1_LE;
%         x_cT_w1_TE = Geo_tier.x_cT_w1_TE;
%         
%         y_cR_w1_LE = Geo_tier.y_cR_w1_LE;
%         z_cR_w1_LE = Geo_tier.z_cR_w1_LE;
%         z_cR_w1_TE = Geo_tier.z_cR_w1_TE;
%         z_cT_w1_LE = Geo_tier.z_cT_w1_LE;
%        
%         %--------------------------- WING 1 ---------------------------------
%         S = S_w1;
%         b = b_w1/2;
%         b_w1_fus_in =(0); % increase in y-direction of span associated to fuselage width (root)
%         b_w1_fus_out = (b_w1_e); % increase in y-direction of span associated to fuselage width (tip)
%         y = [b_w1_fus_in b_w1_fus_out]'/2;
%         c_y = [cR_w1 cT_w1]';
%         le_y_1 = (x_cR_w1_LE - x_cR_w1_TE);
%         le_y_2 = (x_cT_w1_LE - x_cR_w1_TE);
%         le_y = [le_y_1 le_y_2]';
%         diedro = [dihedral_w1 dihedral_w1];
%         NACA_foil = 1;
%         VTP_ms = 0;
%         t_c = 0.18;
%         [x_mesh_w1,y_mesh_w1,z_mesh_w1,CA_w1,MISC_w1] = get_DATA_aero(NACA_foil,S,b,y,c_y,le_y,diedro,VTP_ms,t_c);
%         
%         % Modifies location of Mesh
%         % Wing1
%         X = x_cR_w1_TE;
%         Y = y_offset_w1;
%         Y = 0;
%         center_section = 1; % Flag that determine if symmetry along the x axis needs to be appliesd
%         Z = z_cR_w1_TE;
%         [x_mesh_w1_New, y_mesh_w1_New, z_mesh_w1_New] = get_DATA_New_Location(x_mesh_w1, y_mesh_w1, z_mesh_w1,X,Y,Z,center_section);
%         
%         PLOTTING_UAV.x_mesh_w1_New = x_mesh_w1_New; % WING
%         PLOTTING_UAV.y_mesh_w1_New = y_mesh_w1_New; % WING
%         PLOTTING_UAV.z_mesh_w1_New = z_mesh_w1_New; % WING
% 

        [PLOTTING_UAV,x_w1_LE] = generate_wing_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX);

        % Identifies the locations where to determine information of the fuselaje
        % to ensure that that surfaces take into account the fuselaje
        x_loc = [x_w1_LE];
    case 2 % AC_type = 2 - Conventional 2 surface: wing + HTP + VTP
        
%         l_w1w2 = Geo_tier.l_xac_w1w2;
%         l_w1VTP = Geo_tier.l_xac_w1VTP;
        
%         % Surfaces
%         S_w1 = Geo_tier.S_w1;
%         S_w2 = Geo_tier.S_w2;
%         S_VTP = Geo_tier.S_VTP;
%         
%         % Spans
%         b_w1 = Geo_tier.b_w1;
%         b_w1_e = Geo_tier.b_w1_e;
%         b_w2 = Geo_tier.b_w2;
%         b_w2_e = Geo_tier.b_w2_e;
%         b_VTP = Geo_tier.b_VTP;
%         b_VTP_e = Geo_tier.b_VTP_e;
%         
%         % Root Chords
%         cR_w1 = Geo_tier.cR_w1;
%         cR_w2 = Geo_tier.cR_w2;
%         cR_VTP = Geo_tier.cR_VTP;
%         
%         % Tip Chords
%         cT_w1 = Geo_tier.cT_w1;
%         cT_w2 = Geo_tier.cT_w2;
%         cT_VTP = Geo_tier.cT_VTP;
%         
%         % mean chord
%         c_w1 = Geo_tier.cmac_w1;
%         c_w2 = Geo_tier.cmac_w2;
%         c_VTP = Geo_tier.cmac_VTP;
%         
%         dihedral_w1 = Geo_tier.dihedral_w1;
%         dihedral_w2 = Geo_tier.dihedral_w2;
%         dihedral_VTP = Geo_tier.dihedral_VTP;
%         
%         % Sweep
%         Lambda_LE_w1 = Geo_tier.Lambda_LE_w1;
%         Lambda_TE_w1 = Geo_tier.Lambda_TE_w1;
%         Lambda_c4_w1 = Geo_tier.Lambda_c4_w1;
%         Lambda_c2_w1 = Geo_tier.Lambda_c2_w1;
%         
%         Lambda_LE_w2 = Geo_tier.Lambda_LE_w2;
%         Lambda_TE_w2 = Geo_tier.Lambda_TE_w2;
%         Lambda_c4_w2 = Geo_tier.Lambda_c4_w2;
%         Lambda_c2_w2 = Geo_tier.Lambda_c2_w2;
% 
%         Lambda_LE_VTP = Geo_tier.Lambda_LE_VTP;
%         Lambda_TE_VTP = Geo_tier.Lambda_TE_VTP;
%         Lambda_c4_VTP = Geo_tier.Lambda_c4_VTP;
%         Lambda_c2_VTP = Geo_tier.Lambda_c2_VTP;
% 
%         % %------------------------------GEOMETRY------------------------------%
%         % % posición del borde de ataque de las superficies aerodinámicas
%         x_w1_LE = Geo_tier.x_w1_LE;
%         x_w2_LE = Geo_tier.x_w2_LE;
%         x_VTP_LE = Geo_tier.x_VTP_LE;
%         
%         y_offset_w1 = Geo_tier.y_offset_w1;
%         y_offset_w2 = Geo_tier.y_offset_w2;
%         y_offset_VTP = Geo_tier.y_offset_VTP;
%         
%         x_cR_w1_LE = Geo_tier.x_cR_w1_LE;
%         x_cR_w1_TE = Geo_tier.x_cR_w1_TE;
%         x_cT_w1_LE = Geo_tier.x_cT_w1_LE;
%         x_cT_w1_TE = Geo_tier.x_cT_w1_TE;
%         
%         y_cR_w1_LE = Geo_tier.y_cR_w1_LE;
%         z_cR_w1_LE = Geo_tier.z_cR_w1_LE;
%         z_cR_w1_TE = Geo_tier.z_cR_w1_TE;
%         z_cT_w1_LE = Geo_tier.z_cT_w1_LE;
%         
%         x_cR_w2_LE = Geo_tier.x_cR_w2_LE;
%         x_cR_w2_TE = Geo_tier.x_cR_w2_TE;
%         x_cT_w2_LE = Geo_tier.x_cT_w2_LE;
%         x_cT_w2_TE = Geo_tier.x_cT_w2_TE;
%         
%         y_cR_w2_LE = Geo_tier.y_cR_w2_LE;
%         z_cR_w2_LE = Geo_tier.z_cR_w2_LE;
%         z_cR_w2_TE = Geo_tier.z_cR_w2_TE;
%         z_cT_w2_LE = Geo_tier.z_cT_w2_LE;
% 
%         x_cR_VTP_LE = Geo_tier.x_cR_w2_LE;
%         x_cR_VTP_TE = Geo_tier.x_cR_w2_TE;
%         x_cT_VTP_LE = Geo_tier.x_cT_w2_LE;
%         x_cT_VTP_TE = Geo_tier.x_cT_w2_TE;
%         
%         y_cR_VTP_LE = Geo_tier.y_cR_VTP_LE;
%         z_cR_VTP_LE = Geo_tier.z_cR_VTP_LE;
%         z_cR_VTP_TE = Geo_tier.z_cR_VTP_TE;
%         z_cT_VTP_LE = Geo_tier.z_cT_VTP_LE;
% 
%         %--------------------------- WING 1 ---------------------------------
%         S = S_w1;
%         b = b_w1/2;
%         b_w1_fus_in =(0); % increase in y-direction of span associated to fuselage width (root)
%         b_w1_fus_out = (b_w1_e); % increase in y-direction of span associated to fuselage width (tip)
%         y = [b_w1_fus_in b_w1_fus_out]'/2;
%         c_y = [cR_w1 cT_w1]';
%         le_y_1 = (x_cR_w1_LE - x_cR_w1_TE);
%         le_y_2 = (x_cT_w1_LE - x_cR_w1_TE);
%         le_y = [le_y_1 le_y_2]';
%         diedro = [dihedral_w1 dihedral_w1];
%         NACA_foil = 1;
%         VTP_ms = 0;
%         t_c = 0.18;
%         [x_mesh_w1,y_mesh_w1,z_mesh_w1,CA_w1,MISC_w1] = get_DATA_aero(NACA_foil,S,b,y,c_y,le_y,diedro,VTP_ms,t_c);
%         
%         % Modifies location of Mesh
%         % Wing1
%         X = x_cR_w1_TE;
%         if OUTPUT_read_XLSX.InputGeometry_Data_flags.wing_offset_w1 ==1
%             Y = y_offset_w1;
%         else
%             Y = 0;
%         end
%         center_section = 1; % Flag that determine if symmetry along the x axis needs to be appliesd
%         Z = z_cR_w1_TE;
%         [x_mesh_w1_New, y_mesh_w1_New, z_mesh_w1_New] = get_DATA_New_Location(x_mesh_w1, y_mesh_w1, z_mesh_w1,X,Y,Z,center_section);
%         
%         PLOTTING_UAV.x_mesh_w1_New = x_mesh_w1_New; % WING
%         PLOTTING_UAV.y_mesh_w1_New = y_mesh_w1_New; % WING
%         PLOTTING_UAV.z_mesh_w1_New = z_mesh_w1_New; % WING
%         
%         %--------------------------- WING 2 ---------------------------------
%         S = S_w2;
%         b = b_w2/2;
%         b_w2_fus_in =(0); % increase in y-direction of span associated to fuselage width (root)
%         b_w2_fus_out = (b_w2_e); % increase in y-direction of span associated to fuselage width (tip)
%         y = [b_w2_fus_in b_w2_fus_out]'/2;
%         c_y = [cR_w2 cT_w2]';
%         le_y_1 = (x_cR_w2_LE - x_cR_w2_TE);
%         le_y_2 = (x_cT_w2_LE - x_cR_w2_TE);
%         le_y = [le_y_1 le_y_2]';
%         diedro = [dihedral_w2 dihedral_w2];
%         NACA_foil = 1;
%         VTP_ms = 0;
%         t_c = 0.10  ;
%         [x_mesh_w2, y_mesh_w2, z_mesh_w2,  CA_w2, MISC_w2] = get_DATA_aero(NACA_foil,S,b,y,c_y,le_y,diedro,VTP_ms,t_c);
%         
%         % Wing2
%         X = x_cR_w2_TE;
%         Y = y_offset_w2;
%         center_section = 1;
%         Z = z_cR_w2_TE;
%         [x_mesh_w2_New, y_mesh_w2_New, z_mesh_w2_New] = get_DATA_New_Location(x_mesh_w2, y_mesh_w2, z_mesh_w2,X,Y,Z,center_section);
% 
%         
%         PLOTTING_UAV.x_mesh_w2_New = x_mesh_w2_New; % HTP
%         PLOTTING_UAV.y_mesh_w2_New = y_mesh_w2_New; % HTP
%         PLOTTING_UAV.z_mesh_w2_New = z_mesh_w2_New; % HTP
% 
%         %--------------------------- VTP ---------------------------------
%         S = S_VTP;
%         b = b_VTP;
%         b_VTP_fus_in =(0); % increase in y-direction of span associated to fuselage width (root)
%         b_VTP_fus_out = (b_VTP_e); % increase in y-direction of span associated to fuselage width (tip)
%         y = [b_VTP_fus_in b_VTP_fus_out]';
%         c_y = [cR_VTP cT_VTP]';
%         le_y_1 = (x_cR_VTP_LE - x_cR_VTP_TE);
%         le_y_2 = (x_cT_VTP_LE - x_cT_VTP_TE);
% 
%         le_y = [le_y_1 le_y_2]';
%         dihedral_VTP = 0;
%         diedro = [dihedral_VTP dihedral_VTP];
%         NACA_foil = 1;
%         VTP_ms = 1;
%         t_c = 0.10;
%         [x_mesh_VTP_tmp, y_mesh_VTP_tmp, z_mesh_VTP_tmp,  CA_VTP, MISC_VTP] = get_DATA_aero(NACA_foil,S,b,y,c_y,le_y,diedro,VTP_ms,t_c);
%         
%         x_mesh_VTP =  x_mesh_VTP_tmp;
%         y_mesh_VTP =  y_mesh_VTP_tmp;
%         z_mesh_VTP =  z_mesh_VTP_tmp;
%         
%         % VTP 2
%         X = x_cR_VTP_LE;
%         Y = y_offset_VTP;
%         
%         center_section = 0;
%         Z = z_cR_VTP_TE;
%         [x_mesh_VTP_New_tmp, y_mesh_VTP_New_tmp, z_mesh_VTP_New_tmp] = get_DATA_New_Location(x_mesh_VTP, z_mesh_VTP, y_mesh_VTP,X,Y,Z,center_section);
% 
%         x_mesh_VTP_New =  x_mesh_VTP_New_tmp;
%         y_mesh_VTP_New =  y_mesh_VTP_New_tmp;
%         z_mesh_VTP_New =  z_mesh_VTP_New_tmp;
% 
%         if VTP == 1
%             if twin_VTP == 1
%                 PLOTTING_UAV.x_mesh_VTP1_New = x_mesh_VTP_New; % VTP1
%                 PLOTTING_UAV.y_mesh_VTP1_New = y_mesh_VTP_New; % VTP1
%                 PLOTTING_UAV.z_mesh_VTP1_New = z_mesh_VTP_New; % VTP1
%                 PLOTTING_UAV.x_mesh_VTP2_New = x_mesh_VTP_New; % VTP2
%                 PLOTTING_UAV.y_mesh_VTP2_New = -y_mesh_VTP_New; % VTP2
%                 PLOTTING_UAV.z_mesh_VTP2_New = z_mesh_VTP_New; % VTP2
%         
%             else
%                 PLOTTING_UAV.x_mesh_VTP_New = x_mesh_VTP_New; % VTP
%                 PLOTTING_UAV.y_mesh_VTP_New = y_mesh_VTP_New; % VTP
%                 PLOTTING_UAV.z_mesh_VTP_New = z_mesh_VTP_New; % VTP
%             end
%         end

        % Wing Geometry
        [PLOTTING_UAV,x_w1_LE] = generate_wing_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX);
        % w2 Geometry
        [PLOTTING_UAV, x_w2_LE] = generate_w2_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX);
        % VTP Geometry
        [PLOTTING_UAV, x_VTP_LE] = generate_VTP_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX);

        % Identifies the locations where to determine information of the fuselaje
        % to ensure that that surfaces take into account the fuselaje
        x_loc = [x_w1_LE x_w2_LE x_VTP_LE];

    case 3 % AC_type = 3 - 3 surface: cannard + wing + HTP + VTP

%         l_w1w2 = Geo_tier.l_xac_w1w2;
%         l_w1VTP = Geo_tier.l_xac_w1VTP;
%         
%         % Surfaces
%         S_w1 = Geo_tier.S_w1;
%         S_w2 = Geo_tier.S_w2;
%         S_can = Geo_tier.S_can;
%         S_VTP = Geo_tier.S_VTP;
%         
%         % Spans
%         b_w1 = Geo_tier.b_w1;
%         b_w1_e = Geo_tier.b_w1_e;
%         b_w2 = Geo_tier.b_w2;
%         b_w2_e = Geo_tier.b_w2_e;
%         b_can = Geo_tier.b_can;
%         b_can_e = Geo_tier.b_can_e;
%         b_VTP = Geo_tier.b_VTP;
%         b_VTP_e = Geo_tier.b_VTP_e;
%         
%         % Root Chords
%         cR_w1 = Geo_tier.cR_w1;
%         cR_w2 = Geo_tier.cR_w2;
%         cR_can = Geo_tier.cR_can;
%         cR_VTP = Geo_tier.cR_VTP;
%         
%         % Tip Chords
%         cT_w1 = Geo_tier.cT_w1;
%         cT_w2 = Geo_tier.cT_w2;
%         cT_can = Geo_tier.cT_can;
%         cT_VTP = Geo_tier.cT_VTP;
%         
%         % mean chord
%         c_w1 = Geo_tier.cmac_w1;
%         c_w2 = Geo_tier.cmac_w2;
%         c_can = Geo_tier.cmac_can;
%         c_VTP = Geo_tier.cmac_VTP;
%         
%         dihedral_w1 = Geo_tier.dihedral_w1;
%         dihedral_w2 = Geo_tier.dihedral_w2;
%         dihedral_can = Geo_tier.dihedral_can;
%         dihedral_VTP = Geo_tier.dihedral_VTP;
%         
%         % Sweep
%         Lambda_LE_w1 = Geo_tier.Lambda_LE_w1;
%         Lambda_TE_w1 = Geo_tier.Lambda_TE_w1;
%         Lambda_c4_w1 = Geo_tier.Lambda_c4_w1;
%         Lambda_c2_w1 = Geo_tier.Lambda_c2_w1;
%         
%         Lambda_LE_w2 = Geo_tier.Lambda_LE_w2;
%         Lambda_TE_w2 = Geo_tier.Lambda_TE_w2;
%         Lambda_c4_w2 = Geo_tier.Lambda_c4_w2;
%         Lambda_c2_w2 = Geo_tier.Lambda_c2_w2;
% 
%         Lambda_LE_can = Geo_tier.Lambda_LE_can;
%         Lambda_TE_can = Geo_tier.Lambda_TE_can;
%         Lambda_c4_can = Geo_tier.Lambda_c4_can;
%         Lambda_c2_can = Geo_tier.Lambda_c2_can;
% 
%         Lambda_LE_VTP = Geo_tier.Lambda_LE_VTP;
%         Lambda_TE_VTP = Geo_tier.Lambda_TE_VTP;
%         Lambda_c4_VTP = Geo_tier.Lambda_c4_VTP;
%         Lambda_c2_VTP = Geo_tier.Lambda_c2_VTP;
% 
%         % %------------------------------GEOMETRY------------------------------%
%         % % posición del borde de ataque de las superficies aerodinámicas
%         x_w1_LE = Geo_tier.x_w1_LE;
%         x_w2_LE = Geo_tier.x_w2_LE;
%         x_can_LE = Geo_tier.x_can_LE;
%         x_VTP_LE = Geo_tier.x_VTP_LE;
%         
%         y_offset_w1 = Geo_tier.y_offset_w1;
%         y_offset_w2 = Geo_tier.y_offset_w2;
%         y_offset_can = Geo_tier.y_offset_can;
%         y_offset_VTP = Geo_tier.y_offset_VTP;
%         
%         x_cR_w1_LE = Geo_tier.x_cR_w1_LE;
%         x_cR_w1_TE = Geo_tier.x_cR_w1_TE;
%         x_cT_w1_LE = Geo_tier.x_cT_w1_LE;
%         x_cT_w1_TE = Geo_tier.x_cT_w1_TE;
%         
%         y_cR_w1_LE = Geo_tier.y_cR_w1_LE;
%         z_cR_w1_LE = Geo_tier.z_cR_w1_LE;
%         z_cR_w1_TE = Geo_tier.z_cR_w1_TE;
%         z_cT_w1_LE = Geo_tier.z_cT_w1_LE;
%         
%         x_cR_w2_LE = Geo_tier.x_cR_w2_LE;
%         x_cR_w2_TE = Geo_tier.x_cR_w2_TE;
%         x_cT_w2_LE = Geo_tier.x_cT_w2_LE;
%         x_cT_w2_TE = Geo_tier.x_cT_w2_TE;
%         
%         y_cR_w2_LE = Geo_tier.y_cR_w2_LE;
%         z_cR_w2_LE = Geo_tier.z_cR_w2_LE;
%         z_cR_w2_TE = Geo_tier.z_cR_w2_TE;
%         z_cT_w2_LE = Geo_tier.z_cT_w2_LE;
% 
%         x_cR_can_LE = Geo_tier.x_cR_can_LE;
%         x_cR_can_TE = Geo_tier.x_cR_can_TE;
%         x_cT_can_LE = Geo_tier.x_cT_can_LE;
%         x_cT_can_TE = Geo_tier.x_cT_can_TE;
%         
%         y_cR_can_LE = Geo_tier.y_cR_can_LE;
%         z_cR_can_LE = Geo_tier.z_cR_can_LE;
%         z_cR_can_TE = Geo_tier.z_cR_can_TE;
%         z_cT_can_LE = Geo_tier.z_cT_can_LE;
% 
%         x_cR_VTP_LE = Geo_tier.x_cR_VTP_LE;
%         x_cR_VTP_TE = Geo_tier.x_cR_VTP_TE;
%         x_cT_VTP_LE = Geo_tier.x_cT_VTP_LE;
%         x_cT_VTP_TE = Geo_tier.x_cT_VTP_TE;
%         
%         y_cR_VTP_LE = Geo_tier.y_cR_VTP_LE;
%         z_cR_VTP_LE = Geo_tier.z_cR_VTP_LE;
%         z_cR_VTP_TE = Geo_tier.z_cR_VTP_TE;
%         z_cT_VTP_LE = Geo_tier.z_cT_VTP_LE;
% 
%         %--------------------------- WING 1 ---------------------------------
%         S = S_w1;
%         b = b_w1/2;
%         b_w1_fus_in =(0); % increase in y-direction of span associated to fuselage width (root)
%         b_w1_fus_out = (b_w1_e); % increase in y-direction of span associated to fuselage width (tip)
%         y = [b_w1_fus_in b_w1_fus_out]'/2;
%         c_y = [cR_w1 cT_w1]';
%         le_y_1 = (x_cR_w1_LE - x_cR_w1_TE);
%         le_y_2 = (x_cT_w1_LE - x_cR_w1_TE);
%         le_y = [le_y_1 le_y_2]';
%         diedro = [dihedral_w1 dihedral_w1];
%         NACA_foil = 1;
%         VTP_ms = 0;
%         t_c = 0.18;
%         [x_mesh_w1,y_mesh_w1,z_mesh_w1,CA_w1,MISC_w1] = get_DATA_aero(NACA_foil,S,b,y,c_y,le_y,diedro,VTP_ms,t_c);
%         
%         % Modifies location of Mesh
%         % Wing1
%         X = x_cR_w1_TE;
%         if OUTPUT_read_XLSX.InputGeometry_Data_flags.wing_offset_w1 ==1
%             Y = y_offset_w1;
%         else
%             Y = 0;
%         end
%         center_section = 1; % Flag that determine if symmetry along the x axis needs to be appliesd
%         Z = z_cR_w1_TE;
%         [x_mesh_w1_New, y_mesh_w1_New, z_mesh_w1_New] = get_DATA_New_Location(x_mesh_w1, y_mesh_w1, z_mesh_w1,X,Y,Z,center_section);
%         
%         PLOTTING_UAV.x_mesh_w1_New = x_mesh_w1_New; % WING
%         PLOTTING_UAV.y_mesh_w1_New = y_mesh_w1_New; % WING
%         PLOTTING_UAV.z_mesh_w1_New = z_mesh_w1_New; % WING
% 
%         %--------------------------- WING 2 ---------------------------------
%         S = S_w2;
%         b = b_w2/2;
%         b_w2_fus_in =(0); % increase in y-direction of span associated to fuselage width (root)
%         b_w2_fus_out = (b_w2_e); % increase in y-direction of span associated to fuselage width (tip)
%         y = [b_w2_fus_in b_w2_fus_out]'/2;
%         c_y = [cR_w2 cT_w2]';
%         le_y_1 = (x_cR_w2_LE - x_cR_w2_TE);
%         le_y_2 = (x_cT_w2_LE - x_cR_w2_TE);
%         le_y = [le_y_1 le_y_2]';
%         diedro = [dihedral_w2 dihedral_w2];
%         NACA_foil = 1;
%         VTP_ms = 0;
%         t_c = 0.10  ;
%         [x_mesh_w2, y_mesh_w2, z_mesh_w2,  CA_w2, MISC_w2] = get_DATA_aero(NACA_foil,S,b,y,c_y,le_y,diedro,VTP_ms,t_c);
%         
%         % Wing2
%         X = x_cR_w2_TE;
%         Y = y_offset_w2;
%         center_section = 1;
%         Z = z_cR_w2_TE;
%         [x_mesh_w2_New, y_mesh_w2_New, z_mesh_w2_New] = get_DATA_New_Location(x_mesh_w2, y_mesh_w2, z_mesh_w2,X,Y,Z,center_section);
% 
%         
%         PLOTTING_UAV.x_mesh_w2_New = x_mesh_w2_New; % HTP
%         PLOTTING_UAV.y_mesh_w2_New = y_mesh_w2_New; % HTP
%         PLOTTING_UAV.z_mesh_w2_New = z_mesh_w2_New; % HTP
% 
%         %--------------------------- CANARD ---------------------------------
%         S = S_can;
%         b = b_can/2;
%         b_can_fus_in =(0); % increase in y-direction of span associated to fuselage width (root)
%         b_can_fus_out = (b_can_e); % increase in y-direction of span associated to fuselage width (tip)
%         y = [b_can_fus_in b_can_fus_out]'/2;
%         c_y = [cR_can cT_can]';
%         le_y_1 = (x_cR_can_LE - x_cR_can_TE);
%         le_y_2 = (x_cT_can_LE - x_cR_can_TE);
%         le_y = [le_y_1 le_y_2]';
%         diedro = [dihedral_can dihedral_can];
%         NACA_foil = 1;
%         VTP_ms = 0;
%         t_c = 0.18;
%         [x_mesh_can,y_mesh_can,z_mesh_can,CA_can,MISC_can] = get_DATA_aero(NACA_foil,S,b,y,c_y,le_y,diedro,VTP_ms,t_c);
%         
%         % Modifies location of Mesh
%         % Wing1
%         X = x_cR_can_TE;
%         if OUTPUT_read_XLSX.InputGeometry_Data_flags.wing_offset_can ==1
%             Y = y_offset_can;
%         else
%             Y = 0;
%         end
%         center_section = 1; % Flag that determine if symmetry along the x axis needs to be appliesd
%         Z = z_cR_can_TE;
%         [x_mesh_can_New, y_mesh_can_New, z_mesh_can_New] = get_DATA_New_Location(x_mesh_can, y_mesh_can, z_mesh_can,X,Y,Z,center_section);
%         
%         PLOTTING_UAV.x_mesh_can_New = x_mesh_can_New; % Canard
%         PLOTTING_UAV.y_mesh_can_New = y_mesh_can_New; % Canard
%         PLOTTING_UAV.z_mesh_can_New = z_mesh_can_New; % Canard
% 
%         %--------------------------- VTP ---------------------------------
%         S = S_VTP;
%         b = b_VTP;
%         b_VTP_fus_in =(0); % increase in y-direction of span associated to fuselage width (root)
%         b_VTP_fus_out = (b_VTP_e); % increase in y-direction of span associated to fuselage width (tip)
%         y = [b_VTP_fus_in b_VTP_fus_out]';
%         c_y = [cR_VTP cT_VTP]';
%         le_y_1 = (x_cR_VTP_LE - x_cR_VTP_TE);
%         le_y_2 = (x_cT_VTP_LE - x_cT_VTP_TE);
% 
%         le_y = [le_y_1 le_y_2]';
%         dihedral_VTP = 0;
%         diedro = [dihedral_VTP dihedral_VTP];
%         NACA_foil = 1;
%         VTP_ms = 1;
%         t_c = 0.10;
%         [x_mesh_VTP_tmp, y_mesh_VTP_tmp, z_mesh_VTP_tmp,  CA_VTP, MISC_VTP] = get_DATA_aero(NACA_foil,S,b,y,c_y,le_y,diedro,VTP_ms,t_c);
%         
%         x_mesh_VTP =  x_mesh_VTP_tmp;
%         y_mesh_VTP =  y_mesh_VTP_tmp;
%         z_mesh_VTP =  z_mesh_VTP_tmp;
%         
%         % VTP 2
%         X = x_cR_VTP_LE;
%         Y = y_offset_VTP;
%         
%         center_section = 0;
%         Z = z_cR_VTP_TE;
%         [x_mesh_VTP_New_tmp, y_mesh_VTP_New_tmp, z_mesh_VTP_New_tmp] = get_DATA_New_Location(x_mesh_VTP, z_mesh_VTP, y_mesh_VTP,X,Y,Z,center_section);
% 
%         x_mesh_VTP_New =  x_mesh_VTP_New_tmp;
%         y_mesh_VTP_New =  y_mesh_VTP_New_tmp;
%         z_mesh_VTP_New =  z_mesh_VTP_New_tmp;
% 
%         if VTP == 1
%             if twin_VTP == 1
%                 PLOTTING_UAV.x_mesh_VTP1_New = x_mesh_VTP_New; % VTP1
%                 PLOTTING_UAV.y_mesh_VTP1_New = y_mesh_VTP_New; % VTP1
%                 PLOTTING_UAV.z_mesh_VTP1_New = z_mesh_VTP_New; % VTP1
%                 PLOTTING_UAV.x_mesh_VTP2_New = x_mesh_VTP_New; % VTP2
%                 PLOTTING_UAV.y_mesh_VTP2_New = -y_mesh_VTP_New; % VTP2
%                 PLOTTING_UAV.z_mesh_VTP2_New = z_mesh_VTP_New; % VTP2
%         
%             else
%                 PLOTTING_UAV.x_mesh_VTP_New = x_mesh_VTP_New; % VTP
%                 PLOTTING_UAV.y_mesh_VTP_New = y_mesh_VTP_New; % VTP        % Identifies the locations where to determine information of the fuselaje
%                 PLOTTING_UAV.z_mesh_VTP_New = z_mesh_VTP_New; % VTP
%             end
%         end

        % Wing Geometry
        [PLOTTING_UAV,x_w1_LE] = generate_wing_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX);
        % Canard Geometry
        [PLOTTING_UAV,x_can_LE] = generate_can_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX);
        % w2 Geometry
        [PLOTTING_UAV, x_w2_LE] = generate_w2_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX);
        % VTP Geometry
        [PLOTTING_UAV, x_VTP_LE] = generate_VTP_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX);

        % Identifies the locations where to determine information of the fuselaje
        % to ensure that that surfaces take into account the fuselaje
        x_loc = [x_can_LE x_w1_LE x_w2_LE x_VTP_LE];
        
    case 4 % AC_type = 4 - 2 surface: wing + V-tail

%         cR_w1 = Geo_tier.cR_w1;
%         cR_w2 = Geo_tier.cR_w2;
%         l_w1w2 = Geo_tier.l_xac_w1w2;
%         
%         % Surfaces
%         S_w1 = Geo_tier.S_w1;
%         S_w2 = Geo_tier.S_w2;
%         
%         % Spans
%         b_w1 = Geo_tier.b_w1;
%         b_w1_e = Geo_tier.b_w1_e;
%         b_w2 = Geo_tier.b_w2;
%         b_w2_e = Geo_tier.b_w2_e;
%         
%         % Root Chords
%         cR_w1 = Geo_tier.cR_w1;
%         cR_w2 = Geo_tier.cR_w2;
%         
%         % Tip Chords
%         cT_w1 = Geo_tier.cT_w1;
%         cT_w2 = Geo_tier.cT_w2;
%         % cT_VTP = Geo_tier.cT_v;
%         
%         % mean chord
%         c_w1 = Geo_tier.cmac_w1;
%         c_w2 = Geo_tier.cmac_w2;
%         % c_VTP = Geo_tier.cmac_v;
%         
%         dihedral_w1 = Geo_tier.dihedral_w1;
%         dihedral_w2 = Geo_tier.dihedral_w2;
%         % Geo_tier.dihedral_v = dihedral_v;
%         
%         % Sweep
%         Lambda_LE_w1 = Geo_tier.Lambda_LE_w1;
%         Lambda_TE_w1 = Geo_tier.Lambda_TE_w1;
%         Lambda_c4_w1 = Geo_tier.Lambda_c4_w1;
%         Lambda_c2_w1 = Geo_tier.Lambda_c2_w1;
%         
%         Lambda_LE_w2 = Geo_tier.Lambda_LE_w2;
%         Lambda_TE_w2 = Geo_tier.Lambda_TE_w2;
%         Lambda_c4_w2 = Geo_tier.Lambda_c4_w2;
%         Lambda_c2_w2 = Geo_tier.Lambda_c2_w2;
%         
%         % %------------------------------GEOMETRY------------------------------%
%         % % posición del borde de ataque de las superficies aerodinámicas
%         x_w1_LE = Geo_tier.x_w1_LE;
%         x_w2_LE = Geo_tier.x_w2_LE;
%         
%         y_offset_w1 = Geo_tier.y_offset_w1;
%         y_offset_w2 = Geo_tier.y_offset_w2;
%         
%         x_cR_w1_LE = Geo_tier.x_cR_w1_LE;
%         x_cR_w1_TE = Geo_tier.x_cR_w1_TE;
%         x_cT_w1_LE = Geo_tier.x_cT_w1_LE;
%         x_cT_w1_TE = Geo_tier.x_cT_w1_TE;
%         
%         y_cR_w1_LE = Geo_tier.y_cR_w1_LE;
%         z_cR_w1_LE = Geo_tier.z_cR_w1_LE;
%         z_cR_w1_TE = Geo_tier.z_cR_w1_TE;
%         z_cT_w1_LE = Geo_tier.z_cT_w1_LE;
%         
%         x_cR_w2_LE = Geo_tier.x_cR_w2_LE;
%         x_cR_w2_TE = Geo_tier.x_cR_w2_TE;
%         x_cT_w2_LE = Geo_tier.x_cT_w2_LE;
%         x_cT_w2_TE = Geo_tier.x_cT_w2_TE;
%         
%         y_cR_w2_LE = Geo_tier.y_cR_w2_LE;
%         z_cR_w2_LE = Geo_tier.z_cR_w2_LE;
%         z_cR_w2_TE = Geo_tier.z_cR_w2_TE;
%         z_cT_w2_LE = Geo_tier.z_cT_w2_LE;
%         
%         %--------------------------- WING 1 ---------------------------------
%         S = S_w1;
%         b = b_w1/2;
%         b_w1_fus_in =(0); % increase in y-direction of span associated to fuselage width (root)
%         b_w1_fus_out = (b_w1_e); % increase in y-direction of span associated to fuselage width (tip)
%         y = [b_w1_fus_in b_w1_fus_out]'/2;
%         c_y = [cR_w1 cT_w1]';
%         le_y_1 = (x_cR_w1_LE - x_cR_w1_TE);
%         le_y_2 = (x_cT_w1_LE - x_cR_w1_TE);
% %         le_y_2 = (x_cT_w1_LE - x_cT_w1_TE);
%         
%         le_y = [le_y_1 le_y_2]';
%         diedro = [dihedral_w1 dihedral_w1];
%         NACA_foil = 1;
%         VTP_ms = 0;
%         t_c = 0.18;
%         [x_mesh_w1,y_mesh_w1,z_mesh_w1,CA_w1,MISC_w1] = get_DATA_aero(NACA_foil,S,b,y,c_y,le_y,diedro,VTP_ms,t_c);
%         
%         % Modifies location of Mesh
%         % Wing1
%         X = x_cR_w1_TE;
%         Y = y_offset_w1;
%         center_section = 1; % Flag that determine if symmetry along the x axis needs to be appliesd
%         Z = z_cR_w1_TE;
%         [x_mesh_w1_New, y_mesh_w1_New, z_mesh_w1_New] = get_DATA_New_Location(x_mesh_w1, y_mesh_w1, z_mesh_w1,X,Y,Z,center_section);
%         
%         PLOTTING_UAV.x_mesh_w1_New = x_mesh_w1_New; % WING
%         PLOTTING_UAV.y_mesh_w1_New = y_mesh_w1_New; % WING
%         PLOTTING_UAV.z_mesh_w1_New = z_mesh_w1_New; % WING
%         
%         %--------------------------- WING 2 ---------------------------------
%         S = S_w2;
%         b = b_w2/2;
%         b_w2_fus_in =(0); % increase in y-direction of span associated to fuselage width (root)
%         b_w2_fus_out = (b_w2_e); % increase in y-direction of span associated to fuselage width (tip)
%         y = [b_w2_fus_in b_w2_fus_out]'/2;
%         c_y = [cR_w2 cT_w2]';
%         le_y_1 = (x_cR_w2_LE - x_cR_w2_TE);
%         le_y_2 = (x_cT_w2_LE - x_cR_w2_TE);
%         le_y = [le_y_1 le_y_2]';
%         diedro = [dihedral_w2 dihedral_w2];
%         NACA_foil = 1;
%         VTP_ms = 0;
%         t_c = 0.12  ;
%         [x_mesh_w2, y_mesh_w2, z_mesh_w2,  CA_w2, MISC_w2] = get_DATA_aero(NACA_foil,S,b,y,c_y,le_y,diedro,VTP_ms,t_c);
%         
%         % Wing2
%         X = x_cR_w2_TE;
%         Y = y_offset_w2;
%         center_section = 1;        

        % Wing Geometry
        [PLOTTING_UAV,x_w1_LE] = generate_wing_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX);
        % w2 Geometry
        [PLOTTING_UAV, x_w2_LE] = generate_w2_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX);
        % Identifies the locations where to determine information of the fuselaje
        % to ensure that that surfaces take into account the fuselaje
        x_loc = [x_w1_LE x_w2_LE];
        
    case 5 % AC_type = 5 - 3 surface: cannard + wing + V-tail
        
%         l_w1w2 = Geo_tier.l_xac_w1w2;
%         
%         % Surfaces
%         S_w1 = Geo_tier.S_w1;
%         S_w2 = Geo_tier.S_w2;
%         S_can = Geo_tier.S_can;
%         
%         % Spans
%         b_w1 = Geo_tier.b_w1;
%         b_w1_e = Geo_tier.b_w1_e;
%         b_w2 = Geo_tier.b_w2;
%         b_w2_e = Geo_tier.b_w2_e;
%         b_can = Geo_tier.b_can;
%         b_can_e = Geo_tier.b_can_e;
%         
%         % Root Chords
%         cR_w1 = Geo_tier.cR_w1;
%         cR_w2 = Geo_tier.cR_w2;
%         cR_can = Geo_tier.cR_can;
%         
%         % Tip Chords
%         cT_w1 = Geo_tier.cT_w1;
%         cT_w2 = Geo_tier.cT_w2;
%         cT_can = Geo_tier.cT_can;
%         
%         % mean chord
%         c_w1 = Geo_tier.cmac_w1;
%         c_w2 = Geo_tier.cmac_w2;
%         c_can = Geo_tier.cmac_can;
%         
%         dihedral_w1 = Geo_tier.dihedral_w1;
%         dihedral_w2 = Geo_tier.dihedral_w2;
%         dihedral_can = Geo_tier.dihedral_can;
%         
%         % Sweep
%         Lambda_LE_w1 = Geo_tier.Lambda_LE_w1;
%         Lambda_TE_w1 = Geo_tier.Lambda_TE_w1;
%         Lambda_c4_w1 = Geo_tier.Lambda_c4_w1;
%         Lambda_c2_w1 = Geo_tier.Lambda_c2_w1;
%         
%         Lambda_LE_w2 = Geo_tier.Lambda_LE_w2;
%         Lambda_TE_w2 = Geo_tier.Lambda_TE_w2;
%         Lambda_c4_w2 = Geo_tier.Lambda_c4_w2;
%         Lambda_c2_w2 = Geo_tier.Lambda_c2_w2;
% 
%         Lambda_LE_can = Geo_tier.Lambda_LE_can;
%         Lambda_TE_can = Geo_tier.Lambda_TE_can;
%         Lambda_c4_can = Geo_tier.Lambda_c4_can;
%         Lambda_c2_can = Geo_tier.Lambda_c2_can;
% 
%         % %------------------------------GEOMETRY------------------------------%
%         % % posición del borde de ataque de las superficies aerodinámicas
%         x_w1_LE = Geo_tier.x_w1_LE;
%         x_w2_LE = Geo_tier.x_w2_LE;
%         x_can_LE = Geo_tier.x_can_LE;
%         
%         y_offset_w1 = Geo_tier.y_offset_w1;
%         y_offset_w2 = Geo_tier.y_offset_w2;
%         y_offset_can = Geo_tier.y_offset_can;
%         
%         x_cR_w1_LE = Geo_tier.x_cR_w1_LE;
%         x_cR_w1_TE = Geo_tier.x_cR_w1_TE;
%         x_cT_w1_LE = Geo_tier.x_cT_w1_LE;
%         x_cT_w1_TE = Geo_tier.x_cT_w1_TE;
%         
%         y_cR_w1_LE = Geo_tier.y_cR_w1_LE;
%         z_cR_w1_LE = Geo_tier.z_cR_w1_LE;
%         z_cR_w1_TE = Geo_tier.z_cR_w1_TE;
%         z_cT_w1_LE = Geo_tier.z_cT_w1_LE;
%         
%         x_cR_w2_LE = Geo_tier.x_cR_w2_LE;
%         x_cR_w2_TE = Geo_tier.x_cR_w2_TE;
%         x_cT_w2_LE = Geo_tier.x_cT_w2_LE;
%         x_cT_w2_TE = Geo_tier.x_cT_w2_TE;
%         
%         y_cR_w2_LE = Geo_tier.y_cR_w2_LE;
%         z_cR_w2_LE = Geo_tier.z_cR_w2_LE;
%         z_cR_w2_TE = Geo_tier.z_cR_w2_TE;
%         z_cT_w2_LE = Geo_tier.z_cT_w2_LE;
% 
%         x_cR_can_LE = Geo_tier.x_cR_can_LE;
%         x_cR_can_TE = Geo_tier.x_cR_can_TE;
%         x_cT_can_LE = Geo_tier.x_cT_can_LE;
%         x_cT_can_TE = Geo_tier.x_cT_can_TE;
%         
%         y_cR_can_LE = Geo_tier.y_cR_can_LE;
%         z_cR_can_LE = Geo_tier.z_cR_can_LE;
%         z_cR_can_TE = Geo_tier.z_cR_can_TE;
%         z_cT_can_LE = Geo_tier.z_cT_can_LE;
% 
%         %--------------------------- WING 1 ---------------------------------
%         S = S_w1;
%         b = b_w1/2;
%         b_w1_fus_in =(0); % increase in y-direction of span associated to fuselage width (root)
%         b_w1_fus_out = (b_w1_e); % increase in y-direction of span associated to fuselage width (tip)
%         y = [b_w1_fus_in b_w1_fus_out]'/2;
%         c_y = [cR_w1 cT_w1]';
%         le_y_1 = (x_cR_w1_LE - x_cR_w1_TE);
%         le_y_2 = (x_cT_w1_LE - x_cR_w1_TE);
%         le_y = [le_y_1 le_y_2]';
%         diedro = [dihedral_w1 dihedral_w1];
%         NACA_foil = 1;
%         VTP_ms = 0;
%         t_c = 0.18;
%         [x_mesh_w1,y_mesh_w1,z_mesh_w1,CA_w1,MISC_w1] = get_DATA_aero(NACA_foil,S,b,y,c_y,le_y,diedro,VTP_ms,t_c);
%         
%         % Modifies location of Mesh
%         % Wing1
%         X = x_cR_w1_TE;
%         if OUTPUT_read_XLSX.InputGeometry_Data_flags.wing_offset_w1 ==1
%             Y = y_offset_w1;
%         else
%             Y = 0;
%         end
%         center_section = 1; % Flag that determine if symmetry along the x axis needs to be appliesd
%         Z = z_cR_w1_TE;
%         [x_mesh_w1_New, y_mesh_w1_New, z_mesh_w1_New] = get_DATA_New_Location(x_mesh_w1, y_mesh_w1, z_mesh_w1,X,Y,Z,center_section);
%         
%         PLOTTING_UAV.x_mesh_w1_New = x_mesh_w1_New; % WING
%         PLOTTING_UAV.y_mesh_w1_New = y_mesh_w1_New; % WING
%         PLOTTING_UAV.z_mesh_w1_New = z_mesh_w1_New; % WING
% 
%                 [PLOTTING_UAV,x_loc] = generate_wing_mesh(Geo_tier);
% 
% 
%         %--------------------------- WING 2 ---------------------------------
%         S = S_w2;
%         b = b_w2/2;
%         b_w2_fus_in =(0); % increase in y-direction of span associated to fuselage width (root)
%         b_w2_fus_out = (b_w2_e); % increase in y-direction of span associated to fuselage width (tip)
%         y = [b_w2_fus_in b_w2_fus_out]'/2;
%         c_y = [cR_w2 cT_w2]';
%         le_y_1 = (x_cR_w2_LE - x_cR_w2_TE);
%         le_y_2 = (x_cT_w2_LE - x_cR_w2_TE);
%         le_y = [le_y_1 le_y_2]';
%         diedro = [dihedral_w2 dihedral_w2];
%         NACA_foil = 1;
%         VTP_ms = 0;
%         t_c = 0.10  ;
%         [x_mesh_w2, y_mesh_w2, z_mesh_w2,  CA_w2, MISC_w2] = get_DATA_aero(NACA_foil,S,b,y,c_y,le_y,diedro,VTP_ms,t_c);
%         
%         % Wing2
%         X = x_cR_w2_TE;
%         Y = y_offset_w2;
%         center_section = 1;
%         Z = z_cR_w2_TE;
%         [x_mesh_w2_New, y_mesh_w2_New, z_mesh_w2_New] = get_DATA_New_Location(x_mesh_w2, y_mesh_w2, z_mesh_w2,X,Y,Z,center_section);
% 
%         
%         PLOTTING_UAV.x_mesh_w2_New = x_mesh_w2_New; % HTP
%         PLOTTING_UAV.y_mesh_w2_New = y_mesh_w2_New; % HTP
%         PLOTTING_UAV.z_mesh_w2_New = z_mesh_w2_New; % HTP
% 
%         %--------------------------- CANARD ---------------------------------
%         S = S_can;
%         b = b_can/2;
%         b_can_fus_in =(0); % increase in y-direction of span associated to fuselage width (root)
%         b_can_fus_out = (b_can_e); % increase in y-direction of span associated to fuselage width (tip)
%         y = [b_can_fus_in b_can_fus_out]'/2;
%         c_y = [cR_can cT_can]';
%         le_y_1 = (x_cR_can_LE - x_cR_can_TE);
%         le_y_2 = (x_cT_can_LE - x_cR_can_TE);
%         le_y = [le_y_1 le_y_2]';
%         diedro = [dihedral_can dihedral_can];
%         NACA_foil = 1;
%         VTP_ms = 0;
%         t_c = 0.18;
%         [x_mesh_can,y_mesh_can,z_mesh_can,CA_can,MISC_can] = get_DATA_aero(NACA_foil,S,b,y,c_y,le_y,diedro,VTP_ms,t_c);
%         
%         % Modifies location of Mesh
%         % Wing1
%         X = x_cR_can_TE;
%         if OUTPUT_read_XLSX.InputGeometry_Data_flags.wing_offset_can ==1
%             Y = y_offset_can;
%         else
%             Y = 0;
%         end
%         center_section = 1; % Flag that determine if symmetry along the x axis needs to be appliesd
%         Z = z_cR_can_TE;
%         [x_mesh_can_New, y_mesh_can_New, z_mesh_can_New] = get_DATA_New_Location(x_mesh_can, y_mesh_can, z_mesh_can,X,Y,Z,center_section);
%         
%         PLOTTING_UAV.x_mesh_can_New = x_mesh_can_New; % Canard
%         PLOTTING_UAV.y_mesh_can_New = y_mesh_can_New; % Canard
%         PLOTTING_UAV.z_mesh_can_New = z_mesh_can_New; % Canard

        % Wing Geometry
        [PLOTTING_UAV,x_w1_LE] = generate_wing_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX);
        % Canard Geometry
        [PLOTTING_UAV,x_can_LE] = generate_can_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX);
        % w2 Geometry
        [PLOTTING_UAV, x_w2_LE] = generate_w2_mesh(Geo_tier,PLOTTING_UAV),OUTPUT_read_XLSX;

        % Identifies the locations where to determine information of the fuselaje
        % to ensure that that surfaces take into account the fuselaje
        x_loc = [x_can_LE x_w1_LE x_w2_LE];

    case 6 % AC_type = 6 - 2 surface: cannard + wing + VTP
%         l_w1w2 = Geo_tier.l_xac_w1w2;
%         l_w1VTP = Geo_tier.l_xac_w1VTP;
        
%         % Surfaces
%         S_w1 = Geo_tier.S_w1;
%         S_can = Geo_tier.S_can;
%         S_VTP = Geo_tier.S_VTP;
%         
%         % Spans
%         b_w1 = Geo_tier.b_w1;
%         b_w1_e = Geo_tier.b_w1_e;
%         b_can = Geo_tier.b_can;
%         b_can_e = Geo_tier.b_can_e;
%         b_VTP = Geo_tier.b_VTP;
%         b_VTP_e = Geo_tier.b_VTP_e;
%         
%         % Root Chords
%         cR_w1 = Geo_tier.cR_w1;
%         cR_can = Geo_tier.cR_can;
%         cR_VTP = Geo_tier.cR_VTP;
%         
%         % Tip Chords
%         cT_w1 = Geo_tier.cT_w1;
%         cT_can = Geo_tier.cT_can;
%         cT_VTP = Geo_tier.cT_VTP;
%         
%         % mean chord
%         c_w1 = Geo_tier.cmac_w1;
%         c_can = Geo_tier.cmac_can;
%         c_VTP = Geo_tier.cmac_VTP;
%         
%         dihedral_w1 = Geo_tier.dihedral_w1;
%         dihedral_can = Geo_tier.dihedral_can;
%         dihedral_VTP = Geo_tier.dihedral_VTP;
%         
%         % Sweep
%         Lambda_LE_w1 = Geo_tier.Lambda_LE_w1;
%         Lambda_TE_w1 = Geo_tier.Lambda_TE_w1;
%         Lambda_c4_w1 = Geo_tier.Lambda_c4_w1;
%         Lambda_c2_w1 = Geo_tier.Lambda_c2_w1;
%         
%         Lambda_LE_can = Geo_tier.Lambda_LE_can;
%         Lambda_TE_can = Geo_tier.Lambda_TE_can;
%         Lambda_c4_can = Geo_tier.Lambda_c4_can;
%         Lambda_c2_can = Geo_tier.Lambda_c2_can;
% 
%         Lambda_LE_VTP = Geo_tier.Lambda_LE_VTP;
%         Lambda_TE_VTP = Geo_tier.Lambda_TE_VTP;
%         Lambda_c4_VTP = Geo_tier.Lambda_c4_VTP;
%         Lambda_c2_VTP = Geo_tier.Lambda_c2_VTP;
% 
%         % %------------------------------GEOMETRY------------------------------%
%         % % posición del borde de ataque de las superficies aerodinámicas
%         x_w1_LE = Geo_tier.x_w1_LE;
%         x_can_LE = Geo_tier.x_can_LE;
%         x_VTP_LE = Geo_tier.x_VTP_LE;
%         
%         y_offset_w1 = Geo_tier.y_offset_w1;
%         y_offset_can = Geo_tier.y_offset_can;
%         y_offset_VTP = Geo_tier.y_offset_VTP;
%         
%         x_cR_w1_LE = Geo_tier.x_cR_w1_LE;
%         x_cR_w1_TE = Geo_tier.x_cR_w1_TE;
%         x_cT_w1_LE = Geo_tier.x_cT_w1_LE;
%         x_cT_w1_TE = Geo_tier.x_cT_w1_TE;
%         
%         y_cR_w1_LE = Geo_tier.y_cR_w1_LE;
%         z_cR_w1_LE = Geo_tier.z_cR_w1_LE;
%         z_cR_w1_TE = Geo_tier.z_cR_w1_TE;
%         z_cT_w1_LE = Geo_tier.z_cT_w1_LE;
%         
%         x_cR_can_LE = Geo_tier.x_cR_can_LE;
%         x_cR_can_TE = Geo_tier.x_cR_can_TE;
%         x_cT_can_LE = Geo_tier.x_cT_can_LE;
%         x_cT_can_TE = Geo_tier.x_cT_can_TE;
%         
%         y_cR_can_LE = Geo_tier.y_cR_can_LE;
%         z_cR_can_LE = Geo_tier.z_cR_can_LE;
%         z_cR_can_TE = Geo_tier.z_cR_can_TE;
%         z_cT_can_LE = Geo_tier.z_cT_can_LE;
% 
%         x_cR_VTP_LE = Geo_tier.x_cR_VTP_LE;
%         x_cR_VTP_TE = Geo_tier.x_cR_VTP_TE;
%         x_cT_VTP_LE = Geo_tier.x_cT_VTP_LE;
%         x_cT_VTP_TE = Geo_tier.x_cT_VTP_TE;
%         
%         y_cR_VTP_LE = Geo_tier.y_cR_VTP_LE;
%         z_cR_VTP_LE = Geo_tier.z_cR_VTP_LE;
%         z_cR_VTP_TE = Geo_tier.z_cR_VTP_TE;
%         z_cT_VTP_LE = Geo_tier.z_cT_VTP_LE;
% 
%         %--------------------------- WING 1 ---------------------------------
%         S = S_w1;
%         b = b_w1/2;
%         b_w1_fus_in =(0); % increase in y-direction of span associated to fuselage width (root)
%         b_w1_fus_out = (b_w1_e); % increase in y-direction of span associated to fuselage width (tip)
%         y = [b_w1_fus_in b_w1_fus_out]'/2;
%         c_y = [cR_w1 cT_w1]';
%         le_y_1 = (x_cR_w1_LE - x_cR_w1_TE);
%         le_y_2 = (x_cT_w1_LE - x_cR_w1_TE);
%         le_y = [le_y_1 le_y_2]';
%         diedro = [dihedral_w1 dihedral_w1];
%         NACA_foil = 1;
%         VTP_ms = 0;
%         t_c = 0.18;
%         [x_mesh_w1,y_mesh_w1,z_mesh_w1,CA_w1,MISC_w1] = get_DATA_aero(NACA_foil,S,b,y,c_y,le_y,diedro,VTP_ms,t_c);
%         
%         % Modifies location of Mesh
%         % Wing1
%         X = x_cR_w1_TE;
%         if OUTPUT_read_XLSX.InputGeometry_Data_flags.wing_offset_w1 ==1
%             Y = y_offset_w1;
%         else
%             Y = 0;
%         end
%         center_section = 1; % Flag that determine if symmetry along the x axis needs to be appliesd
%         Z = z_cR_w1_TE;
%         [x_mesh_w1_New, y_mesh_w1_New, z_mesh_w1_New] = get_DATA_New_Location(x_mesh_w1, y_mesh_w1, z_mesh_w1,X,Y,Z,center_section);
%         
%         PLOTTING_UAV.x_mesh_w1_New = x_mesh_w1_New; % WING
%         PLOTTING_UAV.y_mesh_w1_New = y_mesh_w1_New; % WING
%         PLOTTING_UAV.z_mesh_w1_New = z_mesh_w1_New; % WING
% 
%         %--------------------------- CANARD ---------------------------------
%         S = S_can;
%         b = b_can/2;
%         b_can_fus_in =(0); % increase in y-direction of span associated to fuselage width (root)
%         b_can_fus_out = (b_can_e); % increase in y-direction of span associated to fuselage width (tip)
%         y = [b_can_fus_in b_can_fus_out]'/2;
%         c_y = [cR_can cT_can]';
%         le_y_1 = (x_cR_can_LE - x_cR_can_TE);
%         le_y_2 = (x_cT_can_LE - x_cR_can_TE);
%         le_y = [le_y_1 le_y_2]';
%         diedro = [dihedral_can dihedral_can];
%         NACA_foil = 1;
%         VTP_ms = 0;
%         t_c = 0.18;
%         [x_mesh_can,y_mesh_can,z_mesh_can,CA_can,MISC_can] = get_DATA_aero(NACA_foil,S,b,y,c_y,le_y,diedro,VTP_ms,t_c);
%         
%         % Modifies location of Mesh
%         % Wing1
%         X = x_cR_can_TE;
%         if OUTPUT_read_XLSX.InputGeometry_Data_flags.canard_offset_can ==1
%             Y = y_offset_can;
%         else
%             Y = 0;
%         end
%         center_section = 1; % Flag that determine if symmetry along the x axis needs to be appliesd
%         Z = z_cR_can_TE;
%         [x_mesh_can_New, y_mesh_can_New, z_mesh_can_New] = get_DATA_New_Location(x_mesh_can, y_mesh_can, z_mesh_can,X,Y,Z,center_section);
%         
%         PLOTTING_UAV.x_mesh_can_New = x_mesh_can_New; % Canard
%         PLOTTING_UAV.y_mesh_can_New = y_mesh_can_New; % Canard
%         PLOTTING_UAV.z_mesh_can_New = z_mesh_can_New; % Canard
% 
%         %--------------------------- VTP ---------------------------------
%         S = S_VTP;
%         b = b_VTP;
%         b_VTP_fus_in =(0); % increase in y-direction of span associated to fuselage width (root)
%         b_VTP_fus_out = (b_VTP_e); % increase in y-direction of span associated to fuselage width (tip)
%         y = [b_VTP_fus_in b_VTP_fus_out]';
%         c_y = [cR_VTP cT_VTP]';
%         le_y_1 = (x_cR_VTP_LE - x_cR_VTP_TE);
%         le_y_2 = (x_cT_VTP_LE - x_cT_VTP_TE);
% 
%         le_y = [le_y_1 le_y_2]';
%         dihedral_VTP = 0;
%         diedro = [dihedral_VTP dihedral_VTP];
%         NACA_foil = 1;
%         VTP_ms = 1;
%         t_c = 0.10;
%         [x_mesh_VTP_tmp, y_mesh_VTP_tmp, z_mesh_VTP_tmp,  CA_VTP, MISC_VTP] = get_DATA_aero(NACA_foil,S,b,y,c_y,le_y,diedro,VTP_ms,t_c);
%         
%         x_mesh_VTP =  x_mesh_VTP_tmp;
%         y_mesh_VTP =  y_mesh_VTP_tmp;
%         z_mesh_VTP =  z_mesh_VTP_tmp;
%         
%         % VTP 2
%         X = x_cR_VTP_LE;
%         Y = y_offset_VTP;
%         
%         center_section = 0;
%         Z = z_cR_VTP_TE;
%         [x_mesh_VTP_New_tmp, y_mesh_VTP_New_tmp, z_mesh_VTP_New_tmp] = get_DATA_New_Location(x_mesh_VTP, z_mesh_VTP, y_mesh_VTP,X,Y,Z,center_section);
% 
%         x_mesh_VTP_New =  x_mesh_VTP_New_tmp;
%         y_mesh_VTP_New =  y_mesh_VTP_New_tmp;
%         z_mesh_VTP_New =  z_mesh_VTP_New_tmp;
% 
%         if VTP == 1
%             if twin_VTP == 1
%                 PLOTTING_UAV.x_mesh_VTP1_New = x_mesh_VTP_New; % VTP1
%                 PLOTTING_UAV.y_mesh_VTP1_New = y_mesh_VTP_New; % VTP1
%                 PLOTTING_UAV.z_mesh_VTP1_New = z_mesh_VTP_New; % VTP1
%                 PLOTTING_UAV.x_mesh_VTP2_New = x_mesh_VTP_New; % VTP2
%                 PLOTTING_UAV.y_mesh_VTP2_New = -y_mesh_VTP_New; % VTP2
%                 PLOTTING_UAV.z_mesh_VTP2_New = z_mesh_VTP_New; % VTP2
%         
%             else
%                 PLOTTING_UAV.x_mesh_VTP_New = x_mesh_VTP_New; % VTP
%                 PLOTTING_UAV.y_mesh_VTP_New = y_mesh_VTP_New; % VTP
%                 PLOTTING_UAV.z_mesh_VTP_New = z_mesh_VTP_New; % VTP
%             end
%         end

        % Wing Geometry
        [PLOTTING_UAV,x_w1_LE] = generate_wing_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX);
        % Canard Geometry
        [PLOTTING_UAV,x_can_LE] = generate_can_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX);
        % VTP Geometry
        [PLOTTING_UAV, x_VTP_LE] = generate_VTP_mesh(Geo_tier,PLOTTING_UAV,OUTPUT_read_XLSX);

        % Identifies the locations where to determine information of the fuselaje
        % to ensure that that surfaces take into account the fuselaje
        x_loc = [x_can_LE x_w1_LE x_VTP_LE];

end

%--------------------------- FUSELAJE ---------------------------------
% Identifies the locations where to determine information of the fuselaje
% to ensure that that surfaces take into account the fuselaje
% x_loc = [x_w1_LE x_w2_LE];
x_max = Body_Geo.x_max;
y_max = Body_Geo.y_max;
z_max = Body_Geo.z_max;

% Location of aerodynamic Surfaces
for k = 1:length(x_loc)
    Y_vec(k) = interp1(x_max,y_max,x_loc(k),'spline');
    Z_vec(k) = interp1(x_max,z_max,x_loc(k),'spline');
end
PLOTTING_UAV.x_loc = x_loc;
PLOTTING_UAV.Y_vec = Y_vec;
PLOTTING_UAV.Z_vec = Z_vec;

z_max_fus = max(max(meshData{3}));
z_min_fus = min(min(meshData{3}));
PLOTTING_UAV.z_max_fus = z_max_fus;
PLOTTING_UAV.z_min_fus = z_min_fus;

% Propulsive_flags.l_eng = l_eng; %
% Engine Cylinder
eng_dia = OUTPUT_read_XLSX.Propulsive_flags.d_eng; %
eng_length = OUTPUT_read_XLSX.Propulsive_flags.l_eng; %

R_ENG = (eng_dia/2); % radii of disk
H_ENG = eng_length; %heith of disk
SD_ENG = 20; % Side count
[PLOT_ENG] = make_cylinder_special(R_ENG,H_ENG,SD_ENG);

% Engine Cylinder
nac_dia = OUTPUT_read_XLSX.Propulsive_flags.d_nc; %
nac_length = OUTPUT_read_XLSX.Propulsive_flags.l_nc; %

R_NAC = (nac_dia/2); % radii of disk
H_NAC = nac_length; %heith of disk
SD_NAC = 20; % Side count
[PLOT_NAC] = make_cylinder_special(R_NAC,H_NAC,SD_NAC);

%--------------------------- ENGINE ---------------------------------
% Engine Configuration
% Engine location
switch Engine_loc
    case 1 % Engine_loc = 1 - under wings
        % Generates de CYLINDER for the Left and right ENGINES
        % Reference point from the front point
        if n_eng == 2
            X_LOC_ENG = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_eng_ybar1;
            Y_LOC_ENG = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_eng_ybar1;
            Z_LOC_ENG = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_eng_ybar1;
            %         X_LOC_ENG = x_cR_w1_LE + cR_w1/2 - H_ENG/2;
            %         Y_LOC_ENG = (b_w1/2)/2 + R_ENG; % Half wing
            %         Z_LOC_ENG = z_cR_w1_LE + (b_w1_e/2)*tan(dihedral_w1);

            PLOTTING_UAV.Eng_L_x_w1 = PLOT_ENG.PatchData2_X + X_LOC_ENG;
            PLOTTING_UAV.Eng_L_y_w1 = PLOT_ENG.PatchData2_Y - Y_LOC_ENG;
            PLOTTING_UAV.Eng_L_z_w1 = PLOT_ENG.PatchData2_Z + z_cT_w1_LE;
            PLOTTING_UAV.Eng_L_z_w1 = PLOT_ENG.PatchData2_Z + Z_LOC_ENG;

            PLOTTING_UAV.Eng_L_x_w1_s = PLOT_ENG.PatchData1_X + X_LOC_ENG;
            PLOTTING_UAV.Eng_L_y_w1_s = PLOT_ENG.PatchData1_Y - Y_LOC_ENG;
            PLOTTING_UAV.Eng_L_z_w1_s = PLOT_ENG.PatchData1_Z + z_cT_w1_LE;
            PLOTTING_UAV.Eng_L_z_w1_s = PLOT_ENG.PatchData1_Z + Z_LOC_ENG;

            PLOTTING_UAV.Eng_R_x_w1 = PLOT_ENG.PatchData2_X + X_LOC_ENG;
            PLOTTING_UAV.Eng_R_y_w1 = PLOT_ENG.PatchData2_Y + Y_LOC_ENG;
            PLOTTING_UAV.Eng_R_z_w1 = PLOT_ENG.PatchData2_Z + z_cT_w1_LE;
            PLOTTING_UAV.Eng_R_z_w1 = PLOT_ENG.PatchData2_Z + Z_LOC_ENG;

            PLOTTING_UAV.Eng_R_x_w1_s = PLOT_ENG.PatchData1_X + X_LOC_ENG;
            PLOTTING_UAV.Eng_R_y_w1_s = PLOT_ENG.PatchData1_Y + Y_LOC_ENG;
            PLOTTING_UAV.Eng_R_z_w1_s = PLOT_ENG.PatchData1_Z + z_cT_w1_LE;
            PLOTTING_UAV.Eng_R_z_w1_s = PLOT_ENG.PatchData1_Z + Z_LOC_ENG;

            % Nacelle
            X_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_nac_ybar1;
            Y_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_nac_ybar1;
            Z_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_nac_ybar1;

            PLOTTING_UAV.Nac_L_x_w1 = PLOT_NAC.PatchData2_X + X_LOC_NAC;
            PLOTTING_UAV.Nac_L_y_w1 = PLOT_NAC.PatchData2_Y - Y_LOC_NAC;
            PLOTTING_UAV.Nac_L_z_w1 = PLOT_NAC.PatchData2_Z + z_cT_w1_LE;
            PLOTTING_UAV.Nac_L_z_w1 = PLOT_NAC.PatchData2_Z + Z_LOC_NAC;

            PLOTTING_UAV.Nac_L_x_w1_s = PLOT_NAC.PatchData1_X + X_LOC_NAC;
            PLOTTING_UAV.Nac_L_y_w1_s = PLOT_NAC.PatchData1_Y - Y_LOC_NAC;
            PLOTTING_UAV.Nac_L_z_w1_s = PLOT_NAC.PatchData1_Z + z_cT_w1_LE;
            PLOTTING_UAV.Nac_L_z_w1_s = PLOT_NAC.PatchData1_Z + Z_LOC_NAC;

            PLOTTING_UAV.Nac_R_x_w1 = PLOT_NAC.PatchData2_X + X_LOC_NAC;
            PLOTTING_UAV.Nac_R_y_w1 = PLOT_NAC.PatchData2_Y + Y_LOC_NAC;
            PLOTTING_UAV.Nac_R_z_w1 = PLOT_NAC.PatchData2_Z + z_cT_w1_LE;
            PLOTTING_UAV.Nac_R_z_w1 = PLOT_NAC.PatchData2_Z + Z_LOC_NAC;

            PLOTTING_UAV.Nac_R_x_w1_s = PLOT_NAC.PatchData1_X + X_LOC_NAC;
            PLOTTING_UAV.Nac_R_y_w1_s = PLOT_NAC.PatchData1_Y + Y_LOC_NAC;
            PLOTTING_UAV.Nac_R_z_w1_s = PLOT_NAC.PatchData1_Z + z_cT_w1_LE;
            PLOTTING_UAV.Nac_R_z_w1_s = PLOT_NAC.PatchData1_Z + Z_LOC_NAC;
        elseif n_eng == 4

            % First Pair of engines
            X_LOC_ENG = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_eng_ybar1;
            Y_LOC_ENG = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_eng_ybar1;
            Z_LOC_ENG = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_eng_ybar1;
            
            PLOTTING_UAV.Eng_L_x_w1 = PLOT_ENG.PatchData2_X + X_LOC_ENG;
            PLOTTING_UAV.Eng_L_y_w1 = PLOT_ENG.PatchData2_Y - Y_LOC_ENG;
            PLOTTING_UAV.Eng_L_z_w1 = PLOT_ENG.PatchData2_Z + z_cT_w1_LE;
            PLOTTING_UAV.Eng_L_z_w1 = PLOT_ENG.PatchData2_Z + Z_LOC_ENG;

            PLOTTING_UAV.Eng_L_x_w1_s = PLOT_ENG.PatchData1_X + X_LOC_ENG;
            PLOTTING_UAV.Eng_L_y_w1_s = PLOT_ENG.PatchData1_Y - Y_LOC_ENG;
            PLOTTING_UAV.Eng_L_z_w1_s = PLOT_ENG.PatchData1_Z + z_cT_w1_LE;
            PLOTTING_UAV.Eng_L_z_w1_s = PLOT_ENG.PatchData1_Z + Z_LOC_ENG;

            PLOTTING_UAV.Eng_R_x_w1 = PLOT_ENG.PatchData2_X + X_LOC_ENG;
            PLOTTING_UAV.Eng_R_y_w1 = PLOT_ENG.PatchData2_Y + Y_LOC_ENG;
            PLOTTING_UAV.Eng_R_z_w1 = PLOT_ENG.PatchData2_Z + z_cT_w1_LE;
            PLOTTING_UAV.Eng_R_z_w1 = PLOT_ENG.PatchData2_Z + Z_LOC_ENG;

            PLOTTING_UAV.Eng_R_x_w1_s = PLOT_ENG.PatchData1_X + X_LOC_ENG;
            PLOTTING_UAV.Eng_R_y_w1_s = PLOT_ENG.PatchData1_Y + Y_LOC_ENG;
            PLOTTING_UAV.Eng_R_z_w1_s = PLOT_ENG.PatchData1_Z + z_cT_w1_LE;
            PLOTTING_UAV.Eng_R_z_w1_s = PLOT_ENG.PatchData1_Z + Z_LOC_ENG;

            % Second Pair of engines
            X_LOC_ENG = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_eng_ybar2;
            Y_LOC_ENG = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_eng_ybar2;
            Z_LOC_ENG = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_eng_ybar2;
            
            PLOTTING_UAV.Eng_L_x_w2 = PLOT_ENG.PatchData2_X + X_LOC_ENG;
            PLOTTING_UAV.Eng_L_y_w2 = PLOT_ENG.PatchData2_Y - Y_LOC_ENG;
            PLOTTING_UAV.Eng_L_z_w2 = PLOT_ENG.PatchData2_Z + z_cT_w1_LE;
            PLOTTING_UAV.Eng_L_z_w2 = PLOT_ENG.PatchData2_Z + Z_LOC_ENG;

            PLOTTING_UAV.Eng_L_x_w2_s = PLOT_ENG.PatchData1_X + X_LOC_ENG;
            PLOTTING_UAV.Eng_L_y_w2_s = PLOT_ENG.PatchData1_Y - Y_LOC_ENG;
            PLOTTING_UAV.Eng_L_z_w2_s = PLOT_ENG.PatchData1_Z + z_cT_w1_LE;
            PLOTTING_UAV.Eng_L_z_w2_s = PLOT_ENG.PatchData1_Z + Z_LOC_ENG;

            PLOTTING_UAV.Eng_R_x_w2 = PLOT_ENG.PatchData2_X + X_LOC_ENG;
            PLOTTING_UAV.Eng_R_y_w2 = PLOT_ENG.PatchData2_Y + Y_LOC_ENG;
            PLOTTING_UAV.Eng_R_z_w2 = PLOT_ENG.PatchData2_Z + z_cT_w1_LE;
            PLOTTING_UAV.Eng_R_z_w2 = PLOT_ENG.PatchData2_Z + Z_LOC_ENG;

            PLOTTING_UAV.Eng_R_x_w2_s = PLOT_ENG.PatchData1_X + X_LOC_ENG;
            PLOTTING_UAV.Eng_R_y_w2_s = PLOT_ENG.PatchData1_Y + Y_LOC_ENG;
            PLOTTING_UAV.Eng_R_z_w2_s = PLOT_ENG.PatchData1_Z + z_cT_w1_LE;
            PLOTTING_UAV.Eng_R_z_w2_s = PLOT_ENG.PatchData1_Z + Z_LOC_ENG;

            % Nacelle first pair
            X_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_nac_ybar1;
            Y_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_nac_ybar1;
            Z_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_nac_ybar1;

            PLOTTING_UAV.Nac_L_x_w1 = PLOT_NAC.PatchData2_X + X_LOC_NAC;
            PLOTTING_UAV.Nac_L_y_w1 = PLOT_NAC.PatchData2_Y - Y_LOC_NAC;
            PLOTTING_UAV.Nac_L_z_w1 = PLOT_NAC.PatchData2_Z + z_cT_w1_LE;
            PLOTTING_UAV.Nac_L_z_w1 = PLOT_NAC.PatchData2_Z + Z_LOC_NAC;

            PLOTTING_UAV.Nac_L_x_w1_s = PLOT_NAC.PatchData1_X + X_LOC_NAC;
            PLOTTING_UAV.Nac_L_y_w1_s = PLOT_NAC.PatchData1_Y - Y_LOC_NAC;
            PLOTTING_UAV.Nac_L_z_w1_s = PLOT_NAC.PatchData1_Z + z_cT_w1_LE;
            PLOTTING_UAV.Nac_L_z_w1_s = PLOT_NAC.PatchData1_Z + Z_LOC_NAC;

            PLOTTING_UAV.Nac_R_x_w1 = PLOT_NAC.PatchData2_X + X_LOC_NAC;
            PLOTTING_UAV.Nac_R_y_w1 = PLOT_NAC.PatchData2_Y + Y_LOC_NAC;
            PLOTTING_UAV.Nac_R_z_w1 = PLOT_NAC.PatchData2_Z + z_cT_w1_LE;
            PLOTTING_UAV.Nac_R_z_w1 = PLOT_NAC.PatchData2_Z + Z_LOC_NAC;

            PLOTTING_UAV.Nac_R_x_w1_s = PLOT_NAC.PatchData1_X + X_LOC_NAC;
            PLOTTING_UAV.Nac_R_y_w1_s = PLOT_NAC.PatchData1_Y + Y_LOC_NAC;
            PLOTTING_UAV.Nac_R_z_w1_s = PLOT_NAC.PatchData1_Z + z_cT_w1_LE;
            PLOTTING_UAV.Nac_R_z_w1_s = PLOT_NAC.PatchData1_Z + Z_LOC_NAC;

            % Second pair of nacelles
            X_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_nac_ybar2;
            Y_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_nac_ybar2;
            Z_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_nac_ybar2;

            PLOTTING_UAV.Nac_L_x_w2 = PLOT_NAC.PatchData2_X + X_LOC_NAC;
            PLOTTING_UAV.Nac_L_y_w2 = PLOT_NAC.PatchData2_Y - Y_LOC_NAC;
            PLOTTING_UAV.Nac_L_z_w2 = PLOT_NAC.PatchData2_Z + z_cT_w1_LE;
            PLOTTING_UAV.Nac_L_z_w2 = PLOT_NAC.PatchData2_Z + Z_LOC_NAC;

            PLOTTING_UAV.Nac_L_x_w2_s = PLOT_NAC.PatchData1_X + X_LOC_NAC;
            PLOTTING_UAV.Nac_L_y_w2_s = PLOT_NAC.PatchData1_Y - Y_LOC_NAC;
            PLOTTING_UAV.Nac_L_z_w2_s = PLOT_NAC.PatchData1_Z + z_cT_w1_LE;
            PLOTTING_UAV.Nac_L_z_w2_s = PLOT_NAC.PatchData1_Z + Z_LOC_NAC;

            PLOTTING_UAV.Nac_R_x_w2 = PLOT_NAC.PatchData2_X + X_LOC_NAC;
            PLOTTING_UAV.Nac_R_y_w2 = PLOT_NAC.PatchData2_Y + Y_LOC_NAC;
            PLOTTING_UAV.Nac_R_z_w2 = PLOT_NAC.PatchData2_Z + z_cT_w1_LE;
            PLOTTING_UAV.Nac_R_z_w2 = PLOT_NAC.PatchData2_Z + Z_LOC_NAC;

            PLOTTING_UAV.Nac_R_x_w2_s = PLOT_NAC.PatchData1_X + X_LOC_NAC;
            PLOTTING_UAV.Nac_R_y_w2_s = PLOT_NAC.PatchData1_Y + Y_LOC_NAC;
            PLOTTING_UAV.Nac_R_z_w2_s = PLOT_NAC.PatchData1_Z + z_cT_w1_LE;
            PLOTTING_UAV.Nac_R_z_w2_s = PLOT_NAC.PatchData1_Z + Z_LOC_NAC;
        end
    case 2 % Engine_loc = 2 - fuselage front
        if n_eng > 1
            disp('For a front fuselage engine, only one engine can be defined. Please do correct Input Data');
        end
        % Generates de CYLINDER for the Left and right ENGINES
        % Reference point from the front point
        X_LOC_ENG = - H_ENG/2;
        Y_LOC_ENG = 0;
        Z_LOC_ENG = 0;
        
        PLOTTING_UAV.Eng_L_x_w1 = PLOT_ENG.PatchData2_X + X_LOC_ENG;
        PLOTTING_UAV.Eng_L_y_w1 = PLOT_ENG.PatchData2_Y - Y_LOC_ENG;
        PLOTTING_UAV.Eng_L_z_w1 = PLOT_ENG.PatchData2_Z + z_cT_w1_LE;
        PLOTTING_UAV.Eng_L_z_w1 = PLOT_ENG.PatchData2_Z + Z_LOC_ENG;
        
        PLOTTING_UAV.Eng_L_x_w1_s = PLOT_ENG.PatchData1_X + X_LOC_ENG;
        PLOTTING_UAV.Eng_L_y_w1_s = PLOT_ENG.PatchData1_Y - Y_LOC_ENG;
        PLOTTING_UAV.Eng_L_z_w1_s = PLOT_ENG.PatchData1_Z + z_cT_w1_LE;
        PLOTTING_UAV.Eng_L_z_w1_s = PLOT_ENG.PatchData1_Z + Z_LOC_ENG;
        
        % Nacelle
        X_LOC_NAC = - H_ENG/2;
        Y_LOC_NAC = 0;
        Z_LOC_NAC = 0;
        
        PLOTTING_UAV.Nac_L_x_w1 = PLOT_NAC.PatchData2_X + X_LOC_NAC;
        PLOTTING_UAV.Nac_L_y_w1 = PLOT_NAC.PatchData2_Y - Y_LOC_NAC;
        PLOTTING_UAV.Nac_L_z_w1 = PLOT_NAC.PatchData2_Z + z_cT_w1_LE;
        PLOTTING_UAV.Nac_L_z_w1 = PLOT_NAC.PatchData2_Z + Z_LOC_NAC;
        
        PLOTTING_UAV.Nac_L_x_w1_s = PLOT_NAC.PatchData1_X + X_LOC_NAC;
        PLOTTING_UAV.Nac_L_y_w1_s = PLOT_NAC.PatchData1_Y - Y_LOC_NAC;
        PLOTTING_UAV.Nac_L_z_w1_s = PLOT_NAC.PatchData1_Z + z_cT_w1_LE;
        PLOTTING_UAV.Nac_L_z_w1_s = PLOT_NAC.PatchData1_Z + Z_LOC_NAC;
        
    case 3 % Engine_loc = 3 - fuselage rear
        if n_eng > 1
            disp('For a rear fuselage engine, only one engine can be defined. Please do correct Input Data');
        end
        % Generates de CYLINDER for the Left and right ENGINES
        % Reference point from the front point
        X_LOC_ENG = Geo_tier.l_fus - H_ENG/2;
        Y_LOC_ENG = 0;
        Z_LOC_ENG = 0;
        
        PLOTTING_UAV.Eng_L_x_w1 = PLOT_ENG.PatchData2_X + X_LOC_ENG;
        PLOTTING_UAV.Eng_L_y_w1 = PLOT_ENG.PatchData2_Y - Y_LOC_ENG;
        PLOTTING_UAV.Eng_L_z_w1 = PLOT_ENG.PatchData2_Z + z_cT_w1_LE;
        PLOTTING_UAV.Eng_L_z_w1 = PLOT_ENG.PatchData2_Z + Z_LOC_ENG;
        
        PLOTTING_UAV.Eng_L_x_w1_s = PLOT_ENG.PatchData1_X + X_LOC_ENG;
        PLOTTING_UAV.Eng_L_y_w1_s = PLOT_ENG.PatchData1_Y - Y_LOC_ENG;
        PLOTTING_UAV.Eng_L_z_w1_s = PLOT_ENG.PatchData1_Z + z_cT_w1_LE;
        PLOTTING_UAV.Eng_L_z_w1_s = PLOT_ENG.PatchData1_Z + Z_LOC_ENG;

    case 4 % Engine_loc = 4 - wingtips
        % Generates de CYLINDER for the Left and right ENGINES
        % Reference point from the front point
%         X_LOC_ENG = x_cT_w1_LE + cT_w1/4 - H_ENG/2;
%         Y_LOC_ENG = b_w1/2 + R_ENG;
%         Z_LOC_ENG = z_cR_w1_LE + (b_w1_e/2)*tan(dihedral_w1);
        X_LOC_ENG = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_eng_ybar1;
        Y_LOC_ENG = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_eng_ybar1;
        Z_LOC_ENG = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_eng_ybar1;
        
        PLOTTING_UAV.Eng_L_x_w1 = PLOT_ENG.PatchData2_X + X_LOC_ENG;
        PLOTTING_UAV.Eng_L_y_w1 = PLOT_ENG.PatchData2_Y - Y_LOC_ENG;
        PLOTTING_UAV.Eng_L_z_w1 = PLOT_ENG.PatchData2_Z + z_cT_w1_LE;
        PLOTTING_UAV.Eng_L_z_w1 = PLOT_ENG.PatchData2_Z + Z_LOC_ENG;
        
        PLOTTING_UAV.Eng_L_x_w1_s = PLOT_ENG.PatchData1_X + X_LOC_ENG;
        PLOTTING_UAV.Eng_L_y_w1_s = PLOT_ENG.PatchData1_Y - Y_LOC_ENG;
        PLOTTING_UAV.Eng_L_z_w1_s = PLOT_ENG.PatchData1_Z + z_cT_w1_LE;
        PLOTTING_UAV.Eng_L_z_w1_s = PLOT_ENG.PatchData1_Z + Z_LOC_ENG;
        
        PLOTTING_UAV.Eng_R_x_w1 = PLOT_ENG.PatchData2_X + X_LOC_ENG;
        PLOTTING_UAV.Eng_R_y_w1 = PLOT_ENG.PatchData2_Y + Y_LOC_ENG;
        PLOTTING_UAV.Eng_R_z_w1 = PLOT_ENG.PatchData2_Z + z_cT_w1_LE;
        PLOTTING_UAV.Eng_R_z_w1 = PLOT_ENG.PatchData2_Z + Z_LOC_ENG;
        
        PLOTTING_UAV.Eng_R_x_w1_s = PLOT_ENG.PatchData1_X + X_LOC_ENG;
        PLOTTING_UAV.Eng_R_y_w1_s = PLOT_ENG.PatchData1_Y + Y_LOC_ENG;
        PLOTTING_UAV.Eng_R_z_w1_s = PLOT_ENG.PatchData1_Z + z_cT_w1_LE;
        PLOTTING_UAV.Eng_R_z_w1_s = PLOT_ENG.PatchData1_Z + Z_LOC_ENG;

        % Nacelle
        X_LOC_NAC = x_cT_w1_LE + cT_w1/4 - H_ENG/2;
        Y_LOC_NAC = b_w1/2 + R_ENG;
        Z_LOC_NAC = z_cR_w1_LE + (b_w1_e/2)*tan(dihedral_w1);
       % Nacelle
        X_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_nac_ybar1;
        Y_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_nac_ybar1;
        Z_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_nac_ybar1;
        
        PLOTTING_UAV.Nac_L_x_w1 = PLOT_NAC.PatchData2_X + X_LOC_NAC;
        PLOTTING_UAV.Nac_L_y_w1 = PLOT_NAC.PatchData2_Y - Y_LOC_NAC;
        PLOTTING_UAV.Nac_L_z_w1 = PLOT_NAC.PatchData2_Z + z_cT_w1_LE;
        PLOTTING_UAV.Nac_L_z_w1 = PLOT_NAC.PatchData2_Z + Z_LOC_NAC;
        
        PLOTTING_UAV.Nac_L_x_w1_s = PLOT_NAC.PatchData1_X + X_LOC_NAC;
        PLOTTING_UAV.Nac_L_y_w1_s = PLOT_NAC.PatchData1_Y - Y_LOC_NAC;
        PLOTTING_UAV.Nac_L_z_w1_s = PLOT_NAC.PatchData1_Z + z_cT_w1_LE;
        PLOTTING_UAV.Nac_L_z_w1_s = PLOT_NAC.PatchData1_Z + Z_LOC_NAC;
        
        PLOTTING_UAV.Nac_R_x_w1 = PLOT_NAC.PatchData2_X + X_LOC_NAC;
        PLOTTING_UAV.Nac_R_y_w1 = PLOT_NAC.PatchData2_Y + Y_LOC_NAC;
        PLOTTING_UAV.Nac_R_z_w1 = PLOT_NAC.PatchData2_Z + z_cT_w1_LE;
        PLOTTING_UAV.Nac_R_z_w1 = PLOT_NAC.PatchData2_Z + Z_LOC_NAC;
        
        PLOTTING_UAV.Nac_R_x_w1_s = PLOT_NAC.PatchData1_X + X_LOC_NAC;
        PLOTTING_UAV.Nac_R_y_w1_s = PLOT_NAC.PatchData1_Y + Y_LOC_NAC;
        PLOTTING_UAV.Nac_R_z_w1_s = PLOT_NAC.PatchData1_Z + z_cT_w1_LE;
        PLOTTING_UAV.Nac_R_z_w1_s = PLOT_NAC.PatchData1_Z + Z_LOC_NAC;

    case 5 % Engine_loc = Engine_loc = 5 - wingtips for wing and canard configuration n_eng at each side
        % Generates de CYLINDER for the Left and right ENGINES

        % Engine Wing
        X_LOC_ENG = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_eng_ybar1;
        Y_LOC_ENG = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_eng_ybar1;
        Z_LOC_ENG = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_eng_ybar1;
        
        PLOTTING_UAV.Eng_L_x_w1 = PLOT_ENG.PatchData2_X + X_LOC_ENG;
        PLOTTING_UAV.Eng_L_y_w1 = PLOT_ENG.PatchData2_Y - Y_LOC_ENG;
        PLOTTING_UAV.Eng_L_z_w1 = PLOT_ENG.PatchData2_Z + z_cT_w1_LE;
        PLOTTING_UAV.Eng_L_z_w1 = PLOT_ENG.PatchData2_Z + Z_LOC_ENG;
        
        PLOTTING_UAV.Eng_L_x_w1_s = PLOT_ENG.PatchData1_X + X_LOC_ENG;
        PLOTTING_UAV.Eng_L_y_w1_s = PLOT_ENG.PatchData1_Y - Y_LOC_ENG;
        PLOTTING_UAV.Eng_L_z_w1_s = PLOT_ENG.PatchData1_Z + z_cT_w1_LE;
        PLOTTING_UAV.Eng_L_z_w1_s = PLOT_ENG.PatchData1_Z + Z_LOC_ENG;
        
        PLOTTING_UAV.Eng_R_x_w1 = PLOT_ENG.PatchData2_X + X_LOC_ENG;
        PLOTTING_UAV.Eng_R_y_w1 = PLOT_ENG.PatchData2_Y + Y_LOC_ENG;
        PLOTTING_UAV.Eng_R_z_w1 = PLOT_ENG.PatchData2_Z + z_cT_w1_LE;
        PLOTTING_UAV.Eng_R_z_w1 = PLOT_ENG.PatchData2_Z + Z_LOC_ENG;
        
        PLOTTING_UAV.Eng_R_x_w1_s = PLOT_ENG.PatchData1_X + X_LOC_ENG;
        PLOTTING_UAV.Eng_R_y_w1_s = PLOT_ENG.PatchData1_Y + Y_LOC_ENG;
        PLOTTING_UAV.Eng_R_z_w1_s = PLOT_ENG.PatchData1_Z + z_cT_w1_LE;
        PLOTTING_UAV.Eng_R_z_w1_s = PLOT_ENG.PatchData1_Z + Z_LOC_ENG;

        % Engine Canard
        X_LOC_ENG = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_eng_ybar2;
        Y_LOC_ENG = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_eng_ybar2;
        Z_LOC_ENG = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_eng_ybar2;
        
        PLOTTING_UAV.Eng_L_x_w2 = PLOT_ENG.PatchData2_X + X_LOC_ENG;
        PLOTTING_UAV.Eng_L_y_w2 = PLOT_ENG.PatchData2_Y - Y_LOC_ENG;
        PLOTTING_UAV.Eng_L_z_w2 = PLOT_ENG.PatchData2_Z + z_cT_can_LE;
        PLOTTING_UAV.Eng_L_z_w2 = PLOT_ENG.PatchData2_Z + Z_LOC_ENG;
        
        PLOTTING_UAV.Eng_L_x_w2_s = PLOT_ENG.PatchData1_X + X_LOC_ENG;
        PLOTTING_UAV.Eng_L_y_w2_s = PLOT_ENG.PatchData1_Y - Y_LOC_ENG;
        PLOTTING_UAV.Eng_L_z_w2_s = PLOT_ENG.PatchData1_Z + z_cT_can_LE;
        PLOTTING_UAV.Eng_L_z_w2_s = PLOT_ENG.PatchData1_Z + Z_LOC_ENG;
        
        PLOTTING_UAV.Eng_R_x_w2 = PLOT_ENG.PatchData2_X + X_LOC_ENG;
        PLOTTING_UAV.Eng_R_y_w2 = PLOT_ENG.PatchData2_Y + Y_LOC_ENG;
        PLOTTING_UAV.Eng_R_z_w2 = PLOT_ENG.PatchData2_Z + z_cT_can_LE;
        PLOTTING_UAV.Eng_R_z_w2 = PLOT_ENG.PatchData2_Z + Z_LOC_ENG;
        
        PLOTTING_UAV.Eng_R_x_w2_s = PLOT_ENG.PatchData1_X + X_LOC_ENG;
        PLOTTING_UAV.Eng_R_y_w2_s = PLOT_ENG.PatchData1_Y + Y_LOC_ENG;
        PLOTTING_UAV.Eng_R_z_w2_s = PLOT_ENG.PatchData1_Z + z_cT_can_LE;
        PLOTTING_UAV.Eng_R_z_w2_s = PLOT_ENG.PatchData1_Z + Z_LOC_ENG;

        % Nacelle Wing
        X_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_nac_ybar1;
        Y_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_nac_ybar1;
        Z_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_nac_ybar1;
        
        PLOTTING_UAV.Nac_L_x_w1 = PLOT_NAC.PatchData2_X + X_LOC_NAC;
        PLOTTING_UAV.Nac_L_y_w1 = PLOT_NAC.PatchData2_Y - Y_LOC_NAC;
        PLOTTING_UAV.Nac_L_z_w1 = PLOT_NAC.PatchData2_Z + z_cT_w1_LE;
        PLOTTING_UAV.Nac_L_z_w1 = PLOT_NAC.PatchData2_Z + Z_LOC_NAC;
        
        PLOTTING_UAV.Nac_L_x_w1_s = PLOT_NAC.PatchData1_X + X_LOC_NAC;
        PLOTTING_UAV.Nac_L_y_w1_s = PLOT_NAC.PatchData1_Y - Y_LOC_NAC;
        PLOTTING_UAV.Nac_L_z_w1_s = PLOT_NAC.PatchData1_Z + z_cT_w1_LE;
        PLOTTING_UAV.Nac_L_z_w1_s = PLOT_NAC.PatchData1_Z + Z_LOC_NAC;
        
        PLOTTING_UAV.Nac_R_x_w1 = PLOT_NAC.PatchData2_X + X_LOC_NAC;
        PLOTTING_UAV.Nac_R_y_w1 = PLOT_NAC.PatchData2_Y + Y_LOC_NAC;
        PLOTTING_UAV.Nac_R_z_w1 = PLOT_NAC.PatchData2_Z + z_cT_w1_LE;
        PLOTTING_UAV.Nac_R_z_w1 = PLOT_NAC.PatchData2_Z + Z_LOC_NAC;
        
        PLOTTING_UAV.Nac_R_x_w1_s = PLOT_NAC.PatchData1_X + X_LOC_NAC;
        PLOTTING_UAV.Nac_R_y_w1_s = PLOT_NAC.PatchData1_Y + Y_LOC_NAC;
        PLOTTING_UAV.Nac_R_z_w1_s = PLOT_NAC.PatchData1_Z + z_cT_w1_LE;
        PLOTTING_UAV.Nac_R_z_w1_s = PLOT_NAC.PatchData1_Z + Z_LOC_NAC;

        % Nacelle Canard
        X_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_nac_ybar2;
        Y_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_nac_ybar2;
        Z_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_nac_ybar2;
        
        PLOTTING_UAV.Nac_L_x_w2 = PLOT_NAC.PatchData2_X + X_LOC_NAC;
        PLOTTING_UAV.Nac_L_y_w2 = PLOT_NAC.PatchData2_Y - Y_LOC_NAC;
        PLOTTING_UAV.Nac_L_z_w2 = PLOT_NAC.PatchData2_Z + z_cT_can_LE;
        PLOTTING_UAV.Nac_L_z_w2 = PLOT_NAC.PatchData2_Z + Z_LOC_NAC;
        
        PLOTTING_UAV.Nac_L_x_w2_s = PLOT_NAC.PatchData1_X + X_LOC_NAC;
        PLOTTING_UAV.Nac_L_y_w2_s = PLOT_NAC.PatchData1_Y - Y_LOC_NAC;
        PLOTTING_UAV.Nac_L_z_w2_s = PLOT_NAC.PatchData1_Z + z_cT_can_LE;
        PLOTTING_UAV.Nac_L_z_w2_s = PLOT_NAC.PatchData1_Z + Z_LOC_NAC;
        
        PLOTTING_UAV.Nac_R_x_w2 = PLOT_NAC.PatchData2_X + X_LOC_NAC;
        PLOTTING_UAV.Nac_R_y_w2 = PLOT_NAC.PatchData2_Y + Y_LOC_NAC;
        PLOTTING_UAV.Nac_R_z_w2 = PLOT_NAC.PatchData2_Z + z_cT_can_LE;
        PLOTTING_UAV.Nac_R_z_w2 = PLOT_NAC.PatchData2_Z + Z_LOC_NAC;
        
        PLOTTING_UAV.Nac_R_x_w2_s = PLOT_NAC.PatchData1_X + X_LOC_NAC;
        PLOTTING_UAV.Nac_R_y_w2_s = PLOT_NAC.PatchData1_Y + Y_LOC_NAC;
        PLOTTING_UAV.Nac_R_z_w2_s = PLOT_NAC.PatchData1_Z + z_cT_can_LE;
        PLOTTING_UAV.Nac_R_z_w2_s = PLOT_NAC.PatchData1_Z + Z_LOC_NAC;
end


% Engine Configuration
switch Engine_conf
    case 1 % Engine_conf = 1 - pusher prop
        % Prop Disk
        R = Prop_data.D_prop/2; % radii of disk
        H = 0.01; %heith of disk
        SD = 20;
        [PLOT_DISK] = make_cylinder_special(R,H,SD);
        Delta_eng = 0.05;
        
       %--------------------------- Prop ---------------------------------
        % Generates de disk for the Left and right Rotors
        % Reference point from the front point
        X_LOC_PROP = X_LOC_ENG + H_ENG;
        Y_LOC_PROP = Y_LOC_ENG;
        Z_LOC_PROP = Z_LOC_ENG;
        
        % Depending the number of engines
        if n_eng == 1
            PLOTTING_UAV.Prop_L_x_w1 = PLOT_DISK.PatchData2_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_L_y_w1 = PLOT_DISK.PatchData2_Y - Y_LOC_PROP;
            PLOTTING_UAV.Prop_L_z_w1 = PLOT_DISK.PatchData2_Z + Z_LOC_PROP;
            
            PLOTTING_UAV.Prop_L_x_w1_s = PLOT_DISK.PatchData1_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_L_y_w1_s = PLOT_DISK.PatchData1_Y - Y_LOC_PROP;
            PLOTTING_UAV.Prop_L_z_w1_s = PLOT_DISK.PatchData1_Z + Z_LOC_PROP;
            
        elseif n_eng == 2
            PLOTTING_UAV.Prop_L_x_w1 = PLOT_DISK.PatchData2_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_L_y_w1 = PLOT_DISK.PatchData2_Y - Y_LOC_PROP;
            PLOTTING_UAV.Prop_L_z_w1 = PLOT_DISK.PatchData2_Z + Z_LOC_PROP;
            
            PLOTTING_UAV.Prop_L_x_w1_s = PLOT_DISK.PatchData1_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_L_y_w1_s = PLOT_DISK.PatchData1_Y - Y_LOC_PROP;
            PLOTTING_UAV.Prop_L_z_w1_s = PLOT_DISK.PatchData1_Z + Z_LOC_PROP;
            
            PLOTTING_UAV.Prop_R_x_w1 = PLOT_DISK.PatchData2_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_R_y_w1 = PLOT_DISK.PatchData2_Y + Y_LOC_PROP;
            PLOTTING_UAV.Prop_R_z_w1 = PLOT_DISK.PatchData2_Z + Z_LOC_PROP;
            
            PLOTTING_UAV.Prop_R_x_w1_s = PLOT_DISK.PatchData1_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_R_y_w1_s = PLOT_DISK.PatchData1_Y + Y_LOC_PROP;
            PLOTTING_UAV.Prop_R_z_w1_s = PLOT_DISK.PatchData1_Z + Z_LOC_PROP;
       
        elseif n_eng == 4
            % First set of engines
            X_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_nac_ybar1;
            Y_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_nac_ybar1;
            Z_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_nac_ybar1;

            PLOTTING_UAV.Prop_L_x_w1 = PLOT_DISK.PatchData2_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_L_y_w1 = PLOT_DISK.PatchData2_Y - Y_LOC_PROP;
            PLOTTING_UAV.Prop_L_z_w1 = PLOT_DISK.PatchData2_Z + Z_LOC_PROP;

            PLOTTING_UAV.Prop_L_x_w1_s = PLOT_DISK.PatchData1_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_L_y_w1_s = PLOT_DISK.PatchData1_Y - Y_LOC_PROP;
            PLOTTING_UAV.Prop_L_z_w1_s = PLOT_DISK.PatchData1_Z + Z_LOC_PROP;

            PLOTTING_UAV.Prop_R_x_w1 = PLOT_DISK.PatchData2_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_R_y_w1 = PLOT_DISK.PatchData2_Y + Y_LOC_PROP;
            PLOTTING_UAV.Prop_R_z_w1 = PLOT_DISK.PatchData2_Z + Z_LOC_PROP;

            PLOTTING_UAV.Prop_R_x_w1_s = PLOT_DISK.PatchData1_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_R_y_w1_s = PLOT_DISK.PatchData1_Y + Y_LOC_PROP;
            PLOTTING_UAV.Prop_R_z_w1_s = PLOT_DISK.PatchData1_Z + Z_LOC_PROP;

            % Second set of engines
            X_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_nac_ybar2;
            Y_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_nac_ybar2;
            Z_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_nac_ybar2;

            PLOTTING_UAV.Prop_L_x_w2 = PLOT_DISK.PatchData2_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_L_y_w2 = PLOT_DISK.PatchData2_Y - Y_LOC_PROP;
            PLOTTING_UAV.Prop_L_z_w2 = PLOT_DISK.PatchData2_Z + Z_LOC_PROP;

            PLOTTING_UAV.Prop_L_x_w2_s = PLOT_DISK.PatchData1_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_L_y_w2_s = PLOT_DISK.PatchData1_Y - Y_LOC_PROP;
            PLOTTING_UAV.Prop_L_z_w2_s = PLOT_DISK.PatchData1_Z + Z_LOC_PROP;

            PLOTTING_UAV.Prop_R_x_w2 = PLOT_DISK.PatchData2_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_R_y_w2 = PLOT_DISK.PatchData2_Y + Y_LOC_PROP;
            PLOTTING_UAV.Prop_R_z_w2 = PLOT_DISK.PatchData2_Z + Z_LOC_PROP;

            PLOTTING_UAV.Prop_R_x_w2_s = PLOT_DISK.PatchData1_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_R_y_w2_s = PLOT_DISK.PatchData1_Y + Y_LOC_PROP;
            PLOTTING_UAV.Prop_R_z_w2_s = PLOT_DISK.PatchData1_Z + Z_LOC_PROP;
        end

        PLOTTING_UAV.PLOT_DISK = PLOT_DISK; % DATA Disk
        
    case 2 % Engine_conf = 2 - puller prop
        % Prop Disk
        R = Prop_data.D_prop/2; % radii of disk
        H = 0.01; %heith of disk
        SD = 20;
        [PLOT_DISK] = make_cylinder_special(R,H,SD);
        Delta_eng = 0.05;
        
        %--------------------------- Prop ---------------------------------
        % Generates de disk for the Left and right Rotors
        % Reference point from the front point
        X_LOC_PROP = X_LOC_ENG;
        Y_LOC_PROP = Y_LOC_ENG;
        Z_LOC_PROP = Z_LOC_ENG;
        
        % Depending the number of engines
        if n_eng == 1
            PLOTTING_UAV.Prop_L_x_w1 = PLOT_DISK.PatchData2_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_L_y_w1 = PLOT_DISK.PatchData2_Y - Y_LOC_PROP;
            PLOTTING_UAV.Prop_L_z_w1 = PLOT_DISK.PatchData2_Z + Z_LOC_PROP;
            
            PLOTTING_UAV.Prop_L_x_w1_s = PLOT_DISK.PatchData1_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_L_y_w1_s = PLOT_DISK.PatchData1_Y - Y_LOC_PROP;
            PLOTTING_UAV.Prop_L_z_w1_s = PLOT_DISK.PatchData1_Z + Z_LOC_PROP;
            
        elseif n_eng == 2
            PLOTTING_UAV.Prop_L_x_w1 = PLOT_DISK.PatchData2_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_L_y_w1 = PLOT_DISK.PatchData2_Y - Y_LOC_PROP;
            PLOTTING_UAV.Prop_L_z_w1 = PLOT_DISK.PatchData2_Z + Z_LOC_PROP;
            
            PLOTTING_UAV.Prop_L_x_w1_s = PLOT_DISK.PatchData1_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_L_y_w1_s = PLOT_DISK.PatchData1_Y - Y_LOC_PROP;
            PLOTTING_UAV.Prop_L_z_w1_s = PLOT_DISK.PatchData1_Z + Z_LOC_PROP;
            
            PLOTTING_UAV.Prop_R_x_w1 = PLOT_DISK.PatchData2_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_R_y_w1 = PLOT_DISK.PatchData2_Y + Y_LOC_PROP;
            PLOTTING_UAV.Prop_R_z_w1 = PLOT_DISK.PatchData2_Z + Z_LOC_PROP;
            
            PLOTTING_UAV.Prop_R_x_w1_s = PLOT_DISK.PatchData1_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_R_y_w1_s = PLOT_DISK.PatchData1_Y + Y_LOC_PROP;
            PLOTTING_UAV.Prop_R_z_w1_s = PLOT_DISK.PatchData1_Z + Z_LOC_PROP;
        elseif n_eng == 4
            % First set of engines
            X_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_nac_ybar1;
            Y_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_nac_ybar1;
            Z_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_nac_ybar1;

            PLOTTING_UAV.Prop_L_x_w1 = PLOT_DISK.PatchData2_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_L_y_w1 = PLOT_DISK.PatchData2_Y - Y_LOC_PROP;
            PLOTTING_UAV.Prop_L_z_w1 = PLOT_DISK.PatchData2_Z + Z_LOC_PROP;

            PLOTTING_UAV.Prop_L_x_w1_s = PLOT_DISK.PatchData1_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_L_y_w1_s = PLOT_DISK.PatchData1_Y - Y_LOC_PROP;
            PLOTTING_UAV.Prop_L_z_w1_s = PLOT_DISK.PatchData1_Z + Z_LOC_PROP;

            PLOTTING_UAV.Prop_R_x_w1 = PLOT_DISK.PatchData2_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_R_y_w1 = PLOT_DISK.PatchData2_Y + Y_LOC_PROP;
            PLOTTING_UAV.Prop_R_z_w1 = PLOT_DISK.PatchData2_Z + Z_LOC_PROP;

            PLOTTING_UAV.Prop_R_x_w1_s = PLOT_DISK.PatchData1_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_R_y_w1_s = PLOT_DISK.PatchData1_Y + Y_LOC_PROP;
            PLOTTING_UAV.Prop_R_z_w1_s = PLOT_DISK.PatchData1_Z + Z_LOC_PROP;

            % Second set of engines
            X_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_nac_ybar2;
            Y_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_nac_ybar2;
            Z_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_nac_ybar2;

            PLOTTING_UAV.Prop_L_x_w2 = PLOT_DISK.PatchData2_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_L_y_w2 = PLOT_DISK.PatchData2_Y - Y_LOC_PROP;
            PLOTTING_UAV.Prop_L_z_w2 = PLOT_DISK.PatchData2_Z + Z_LOC_PROP;

            PLOTTING_UAV.Prop_L_x_w2_s = PLOT_DISK.PatchData1_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_L_y_w2_s = PLOT_DISK.PatchData1_Y - Y_LOC_PROP;
            PLOTTING_UAV.Prop_L_z_w2_s = PLOT_DISK.PatchData1_Z + Z_LOC_PROP;

            PLOTTING_UAV.Prop_R_x_w2 = PLOT_DISK.PatchData2_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_R_y_w2 = PLOT_DISK.PatchData2_Y + Y_LOC_PROP;
            PLOTTING_UAV.Prop_R_z_w2 = PLOT_DISK.PatchData2_Z + Z_LOC_PROP;

            PLOTTING_UAV.Prop_R_x_w2_s = PLOT_DISK.PatchData1_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_R_y_w2_s = PLOT_DISK.PatchData1_Y + Y_LOC_PROP;
            PLOTTING_UAV.Prop_R_z_w2_s = PLOT_DISK.PatchData1_Z + Z_LOC_PROP;
        end
        
        PLOTTING_UAV.PLOT_DISK = PLOT_DISK; % DATA Disk
        
    case 3 % Engine_conf = 3 - turbine
        
        % Prop Disk
        R = Prop_data.D_prop/2; % radii of disk
        H = 0.01; %heith of disk
        SD = 1;
        [PLOT_DISK] = make_cylinder_special(R,H,SD);
        Delta_eng = 0.05;
        
        %--------------------------- Prop ---------------------------------
        % Generates de disk for the Left and right Rotors
        % Reference point from the front point
        X_LOC_PROP = X_LOC_ENG;
        Y_LOC_PROP = Y_LOC_ENG;
        Z_LOC_PROP = Z_LOC_ENG;
        
        % Depending the number of engines
        if n_eng == 1
            PLOTTING_UAV.Prop_L_x_w1 = PLOT_DISK.PatchData2_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_L_y_w1 = PLOT_DISK.PatchData2_Y - Y_LOC_PROP;
            PLOTTING_UAV.Prop_L_z_w1 = PLOT_DISK.PatchData2_Z + Z_LOC_PROP;
            
            PLOTTING_UAV.Prop_L_x_w1_s = PLOT_DISK.PatchData1_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_L_y_w1_s = PLOT_DISK.PatchData1_Y - Y_LOC_PROP;
            PLOTTING_UAV.Prop_L_z_w1_s = PLOT_DISK.PatchData1_Z + Z_LOC_PROP;
            
        elseif n_eng == 2
            PLOTTING_UAV.Prop_L_x_w1 = PLOT_DISK.PatchData2_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_L_y_w1 = PLOT_DISK.PatchData2_Y - Y_LOC_PROP;
            PLOTTING_UAV.Prop_L_z_w1 = PLOT_DISK.PatchData2_Z + Z_LOC_PROP;
            
            PLOTTING_UAV.Prop_L_x_w1_s = PLOT_DISK.PatchData1_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_L_y_w1_s = PLOT_DISK.PatchData1_Y - Y_LOC_PROP;
            PLOTTING_UAV.Prop_L_z_w1_s = PLOT_DISK.PatchData1_Z + Z_LOC_PROP;
            
            PLOTTING_UAV.Prop_R_x_w1 = PLOT_DISK.PatchData2_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_R_y_w1 = PLOT_DISK.PatchData2_Y + Y_LOC_PROP;
            PLOTTING_UAV.Prop_R_z_w1 = PLOT_DISK.PatchData2_Z + Z_LOC_PROP;
            
            PLOTTING_UAV.Prop_R_x_w1_s = PLOT_DISK.PatchData1_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_R_y_w1_s = PLOT_DISK.PatchData1_Y + Y_LOC_PROP;
            PLOTTING_UAV.Prop_R_z_w1_s = PLOT_DISK.PatchData1_Z + Z_LOC_PROP;
        elseif n_eng == 4
            % First set of engines
            X_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_nac_ybar1;
            Y_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_nac_ybar1;
            Z_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_nac_ybar1;

            PLOTTING_UAV.Prop_L_x_w1 = PLOT_DISK.PatchData2_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_L_y_w1 = PLOT_DISK.PatchData2_Y - Y_LOC_PROP;
            PLOTTING_UAV.Prop_L_z_w1 = PLOT_DISK.PatchData2_Z + Z_LOC_PROP;

            PLOTTING_UAV.Prop_L_x_w1_s = PLOT_DISK.PatchData1_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_L_y_w1_s = PLOT_DISK.PatchData1_Y - Y_LOC_PROP;
            PLOTTING_UAV.Prop_L_z_w1_s = PLOT_DISK.PatchData1_Z + Z_LOC_PROP;

            PLOTTING_UAV.Prop_R_x_w1 = PLOT_DISK.PatchData2_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_R_y_w1 = PLOT_DISK.PatchData2_Y + Y_LOC_PROP;
            PLOTTING_UAV.Prop_R_z_w1 = PLOT_DISK.PatchData2_Z + Z_LOC_PROP;

            PLOTTING_UAV.Prop_R_x_w1_s = PLOT_DISK.PatchData1_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_R_y_w1_s = PLOT_DISK.PatchData1_Y + Y_LOC_PROP;
            PLOTTING_UAV.Prop_R_z_w1_s = PLOT_DISK.PatchData1_Z + Z_LOC_PROP;

            % Second set of engines
            X_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.x_nac_ybar2;
            Y_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.y_nac_ybar2;
            Z_LOC_NAC = OUTPUT_read_XLSX.InputGeometry_Data_flags.z_nac_ybar2;

            PLOTTING_UAV.Prop_L_x_w2 = PLOT_DISK.PatchData2_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_L_y_w2 = PLOT_DISK.PatchData2_Y - Y_LOC_PROP;
            PLOTTING_UAV.Prop_L_z_w2 = PLOT_DISK.PatchData2_Z + Z_LOC_PROP;

            PLOTTING_UAV.Prop_L_x_w2_s = PLOT_DISK.PatchData1_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_L_y_w2_s = PLOT_DISK.PatchData1_Y - Y_LOC_PROP;
            PLOTTING_UAV.Prop_L_z_w2_s = PLOT_DISK.PatchData1_Z + Z_LOC_PROP;

            PLOTTING_UAV.Prop_R_x_w2 = PLOT_DISK.PatchData2_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_R_y_w2 = PLOT_DISK.PatchData2_Y + Y_LOC_PROP;
            PLOTTING_UAV.Prop_R_z_w2 = PLOT_DISK.PatchData2_Z + Z_LOC_PROP;

            PLOTTING_UAV.Prop_R_x_w2_s = PLOT_DISK.PatchData1_X + X_LOC_PROP;
            PLOTTING_UAV.Prop_R_y_w2_s = PLOT_DISK.PatchData1_Y + Y_LOC_PROP;
            PLOTTING_UAV.Prop_R_z_w2_s = PLOT_DISK.PatchData1_Z + Z_LOC_PROP;
        end
      PLOTTING_UAV.PLOT_DISK = PLOT_DISK; % DATA Disk        
end

PLOTTING_UAV.meshData = meshData; % FUSELAJE
PLOTTING_UAV.PLOT_ENG = PLOT_ENG; % DATA ENGINE
