% Function that generates the geometry of the different elements
function [MESH_DATA,Fig] = Generation_Geometry_AC(Geo_tier,Plot_Options,Body_Geo,meshData_fus,Prop_data,Fig,AC_CONFIGURATION,COLOR_scheme,case_AC,Performance,CASE_fuse,ESCALADO,SF)

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

[XFLR5_DATA,XFLR5_file,STL_file] = Generation_data_fuselage_EMERGENTIA;

[Body_Geo,meshData_fus] = Generation_Fuselage_Data_old(Geo_tier,XFLR5_DATA,CASE_fuse,ESCALADO,XFLR5_file,STL_file,SF); % Defines Propulsion DATA


%% Defines colors RGB
color_fus = COLOR_scheme.color_fus;
color_w1 = COLOR_scheme.color_w1;
color_w2 = COLOR_scheme.color_w2;
color_prop = COLOR_scheme.color_prop;
color_vtp = COLOR_scheme.color_vtp;
color_eng = COLOR_scheme.color_eng;

if W1 ==1
    %--------------------------- WING 1 ---------------------------------
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
    
    % %------------------------------GEOMETRY------------------------------%
    % % posición del borde de ataque de las superficies aerodinámicas
    x_w1_LE = Geo_tier.x_w1_LE;
    y_offset_w1 = Geo_tier.y_offset_w1;
    x_cR_w1_LE = Geo_tier.x_cR_w1_LE;
    x_cR_w1_TE = Geo_tier.x_cR_w1_TE;
    x_cT_w1_LE = Geo_tier.x_cT_w1_LE;
    x_cT_w1_TE = Geo_tier.x_cT_w1_TE;
    y_cR_w1_LE = Geo_tier.y_cR_w1_LE;
    z_cR_w1_LE = Geo_tier.z_cR_w1_LE;
    z_cR_w1_TE = Geo_tier.z_cR_w1_TE;
    z_cT_w1_LE = Geo_tier.z_cT_w1_LE;
    % Data for obtaining aerodynamic surface
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
    % Changes lovation of Mesh to real relative location with respect to
    % the real design
    % Wing1
    X = x_cR_w1_TE;
    Y = y_offset_w1;
    center_section = 1; % Flag that determine if symmetry along the x axis needs to be appliesd
    Z = z_cR_w1_TE;
    [x_mesh_w1_New, y_mesh_w1_New, z_mesh_w1_New] = get_DATA_New_Location(x_mesh_w1, y_mesh_w1, z_mesh_w1,X,Y,Z,center_section);
    
    % Asigns the color
    C_w1(:,:,1) = color_w1(1)*ones(size(x_mesh_w1_New));
    C_w1(:,:,2) = color_w1(2)*ones(size(y_mesh_w1_New));
    C_w1(:,:,3) = color_w1(3)*ones(size(z_mesh_w1_New));
        
    % Stores the meshing Data
    MESH_DATA.stored_mesh_w1.x_mesh_w1 = x_mesh_w1_New;
    MESH_DATA.stored_mesh_w1.y_mesh_w1 = y_mesh_w1_New;
    MESH_DATA.stored_mesh_w1.z_mesh_w1 = z_mesh_w1_New;
    MESH_DATA.stored_mesh_w1.C_w1 = C_w1;
    
end

if HTP ==1
    %--------------------------- HTP ---------------------------------
    S_w2 = Geo_tier.S_w2;
    % Spans
    b_w2 = Geo_tier.b_w2;
    b_w2_e = Geo_tier.b_w2_e;
    % Root Chords
    cR_w2 = Geo_tier.cR_w2;
    % Tip Chords
    cT_w2 = Geo_tier.cT_w2;
    % mean chord
    c_w2 = Geo_tier.cmac_w2;
    dihedral_w2 = Geo_tier.dihedral_w2;
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
    % Data for obtaining aerodynamic surface
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
    t_c = 0.18;
    [x_mesh_w2,y_mesh_w2,z_mesh_w2,CA_w2,MISC_w2] = get_DATA_aero(NACA_foil,S,b,y,c_y,le_y,diedro,VTP_ms,t_c);
    % Changes lovation of Mesh to real relative location with respect to
    % the real design
    % Wing1
    X = x_cR_w2_TE;
    Y = y_offset_w2;
    center_section = 1; % Flag that determine if symmetry along the x axis needs to be appliesd
    Z = z_cR_w2_TE;
    [x_mesh_w2_New, y_mesh_w2_New, z_mesh_w2_New] = get_DATA_New_Location(x_mesh_w2, y_mesh_w2, z_mesh_w2,X,Y,Z,center_section);
       
    % Asigns the color
    C_w2(:,:,1) = color_w2(1)*ones(size(x_mesh_w2_New));
    C_w2(:,:,2) = color_w2(2)*ones(size(y_mesh_w2_New));
    C_w2(:,:,3) = color_w2(3)*ones(size(z_mesh_w2_New));
        
    % Stores the meshing Data
    MESH_DATA.stored_mesh_w1.x_mesh_w2 = x_mesh_w2_New;
    MESH_DATA.stored_mesh_w1.y_mesh_w2 = y_mesh_w2_New;
    MESH_DATA.stored_mesh_w1.z_mesh_w2 = z_mesh_w2_New;
    MESH_DATA.stored_mesh_w1.C_w2 = C_w2;
    
end

if Vee ==1
    %--------------------------- Vee-tail ---------------------------------
    S_w2 = Geo_tier.S_w2;
    % Spans
    b_w2 = Geo_tier.b_w2;
    b_w2_e = Geo_tier.b_w2_e;
    % Root Chords
    cR_w2 = Geo_tier.cR_w2;
    % Tip Chords
    cT_w2 = Geo_tier.cT_w2;
    % mean chord
    c_w2 = Geo_tier.cmac_w2;
    dihedral_w2 = Geo_tier.dihedral_w2;
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
    % Data for obtaining aerodynamic surface
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
    t_c = 0.18;
    [x_mesh_w2,y_mesh_w2,z_mesh_w2,CA_w2,MISC_w2] = get_DATA_aero(NACA_foil,S,b,y,c_y,le_y,diedro,VTP_ms,t_c);
    % Changes lovation of Mesh to real relative location with respect to
    % the real design
    % Wing1
    X = x_cR_w2_TE;
    Y = y_offset_w2;
    center_section = 1; % Flag that determine if symmetry along the x axis needs to be appliesd
    Z = z_cR_w2_TE;
    [x_mesh_w2_New, y_mesh_w2_New, z_mesh_w2_New] = get_DATA_New_Location(x_mesh_w2, y_mesh_w2, z_mesh_w2,X,Y,Z,center_section);
       
    % Asigns the color
    C_w2(:,:,1) = color_w2(1)*ones(size(x_mesh_w2_New));
    C_w2(:,:,2) = color_w2(2)*ones(size(y_mesh_w2_New));
    C_w2(:,:,3) = color_w2(3)*ones(size(z_mesh_w2_New));
        
    % Stores the meshing Data
    MESH_DATA.stored_mesh_vee.x_mesh_w2 = x_mesh_vee_New;
    MESH_DATA.stored_mesh_vee.y_mesh_w2 = y_mesh_vee_New;
    MESH_DATA.stored_mesh_vee.z_mesh_w2 = z_mesh_vee_New;
    MESH_DATA.stored_mesh_vee.C_w2 = C_vee;
    
end

if VTP ==1
    %--------------------------- VTP ---------------------------------
    S_VTP = Geo_tier.S_VTP;
    % Spans
    b_VTP = Geo_tier.b_VTP;
    b_VTP_e = Geo_tier.b_VTP_e;
    % Root Chords
    cR_VTP = Geo_tier.cR_VTP;
    % Tip Chords
    cT_w2 = Geo_tier.cT_w2;
    % mean chord
    cT_VTP = Geo_tier.cT_VTP;
    dihedral_VTP = Geo_tier.dihedral_VTP;
    % Sweep
    Lambda_LE_VTP = Geo_tier.Lambda_LE_VTP;
    Lambda_TE_VTP = Geo_tier.Lambda_TE_VTP;
    Lambda_c4_VTP = Geo_tier.Lambda_c4_VTP;
    Lambda_c2_VTP = Geo_tier.Lambda_c2_VTP;
    
    %------------------------------GEOMETRY------------------------------%
    % posición del borde de ataque de las superficies aerodinámicas
    x_VTP_LE = Geo_tier.x_VTP_LE;
    y_offset_VTP = Geo_tier.y_offset_VTP;
    x_cR_VTP_LE = Geo_tier.x_cR_VTP_LE;
    x_cR_VTP_TE = Geo_tier.x_cR_VTP_TE;
    x_cT_VTP_LE = Geo_tier.x_cT_VTP_LE;
    x_cT_VTP_TE = Geo_tier.x_cT_VTP_TE;
    y_cR_VTP_LE = Geo_tier.y_cR_VTP_LE;
    z_cR_VTP_LE = Geo_tier.z_cR_VTP_LE;
    z_cR_VTP_TE = Geo_tier.z_cR_VTP_TE;
    z_cT_VTP_LE = Geo_tier.z_cT_VTP_LE;
    % Data for obtaining aerodynamic surface
    S = S_VTP;
    b = b_VTP;
    b_VTP_fus_in =(0); % increase in y-direction of span associated to fuselage width (root)
    b_VTP_fus_out = (b_VTP_e); % increase in y-direction of span associated to fuselage width (tip)
    y = [b_VTP_fus_in b_VTP_fus_out]';
    c_y = [cR_VTP cT_VTP]';
    le_y_1 = (x_cR_VTP_LE - x_cR_VTP_TE);
    le_y_2 = (x_cT_VTP_LE - x_cR_VTP_TE);
    le_y = [le_y_1 le_y_2]';
    diedro = [dihedral_VTP dihedral_VTP];
    NACA_foil = 1;
    VTP_ms = 1;
    t_c = 0.18;
    [x_mesh_VTP,y_mesh_VTP,z_mesh_VTP,CA_VTP,MISC_VTP] = get_DATA_aero(NACA_foil,S,b,y,c_y,le_y,diedro,VTP_ms,t_c);
    % Changes lovation of Mesh to real relative location with respect to
    % the real design
    % Wing1
    X = x_cR_VTP_TE;
    Y = y_offset_VTP;
    center_section = 0; % Flag that determine if symmetry along the x axis needs to be appliesd
    Z = z_cR_VTP_TE;
    [x_mesh_VTP_New, z_mesh_VTP_New, y_mesh_VTP_New] = get_DATA_New_Location(x_mesh_VTP, y_mesh_VTP, z_mesh_VTP,X,Y,Z,center_section);
       
    % Asigns the color
    C_VTP(:,:,1) = color_vtp(1)*ones(size(x_mesh_VTP_New));
    C_VTP(:,:,2) = color_vtp(2)*ones(size(y_mesh_VTP_New));
    C_VTP(:,:,3) = color_vtp(3)*ones(size(z_mesh_VTP_New));
        
    if VTP == 1
        if twin_VTP == 1
            x_mesh_VTP1_New = x_mesh_VTP_New;
            y_mesh_VTP1_New = y_mesh_VTP_New;
            z_mesh_VTP1_New = z_mesh_VTP_New;
            MESH_DATA.stored_mesh_VTP.x_mesh_VTP1_New = x_mesh_VTP1_New; % VTP1
            MESH_DATA.stored_mesh_VTP.y_mesh_VTP1_New = y_mesh_VTP1_New; % VTP1
            MESH_DATA.stored_mesh_VTP.z_mesh_VTP1_New = z_mesh_VTP1_New; % VTP1
            
            x_mesh_VTP2_New = x_mesh_VTP_New;
            y_mesh_VTP2_New = y_mesh_VTP_New;
            z_mesh_VTP2_New = z_mesh_VTP_New;          
            MESH_DATA.stored_mesh_VTP.x_mesh_VTP2_New = x_mesh_VTP2_New; % VTP2
            MESH_DATA.stored_mesh_VTP.y_mesh_VTP2_New = y_mesh_VTP2_New; % VTP2
            MESH_DATA.stored_mesh_VTP.z_mesh_VTP2_New = z_mesh_VTP2_New; % VTP2
            MESH_DATA.stored_mesh_VTP.C_w2 = C_VTP;
        else
            MESH_DATA.stored_mesh_VTP.x_mesh_VTP_New = x_mesh_VTP_New; % VTP
            MESH_DATA.stored_mesh_VTP.y_mesh_VTP_New = y_mesh_VTP_New; % VTP
            MESH_DATA.stored_mesh_VTP.z_mesh_VTP_New = z_mesh_VTP_New; % VTP
            MESH_DATA.stored_mesh_VTP.C_w2 = C_VTP;
        end
    end
           
    % Stores the meshing Data
    MESH_DATA.stored_mesh_VTP.x_mesh_w2 = x_mesh_VTP_New;
    MESH_DATA.stored_mesh_VTP.y_mesh_w2 = y_mesh_VTP_New;
    MESH_DATA.stored_mesh_VTP.z_mesh_w2 = z_mesh_VTP_New;
    MESH_DATA.stored_mesh_VTP.C_w2 = C_VTP;    
end

testing = 1;
if testing == 1
    % plotting example
    Fig = Fig +1;
    figure(Fig)
    sol_mesh_w1 = mesh(x_mesh_w1_New,y_mesh_w1_New,z_mesh_w1_New,C_w1)
    hold on
    pause
    sol_mesh_w2 = mesh(x_mesh_w2_New,y_mesh_w2_New,z_mesh_w2_New,C_w2)
    pause
    sol_mesh_VTP = mesh(x_mesh_VTP1_New,y_mesh_VTP1_New,z_mesh_VTP1_New,C_VTP)
    hold off
end
