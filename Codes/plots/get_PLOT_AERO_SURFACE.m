function [x_mesh,y_mesh,z_mesh,CA,MISC] = get_PLOT_AERO_SURFACE(S,b,cR_w,cT_w,x_w1R_LE,x_w1R_TE,x_w1T_LE,x_w1T_TE,NACA_foil,VTP,dihedral_w)

%--------------------------- WING 1 ---------------------------------
S = S_w1;
b = b_w1/2;
% b_w1_fus_in =(0 + Y_vec(1)); % increase in y-direction of span associated to fuselage width (root)
% b_w1_fus_out = (b_w1 + Y_vec(1)); % increase in y-direction of span associated to fuselage width (tip)
b_w1_fus_in =(0); % increase in y-direction of span associated to fuselage width (root)
b_w1_fus_out = (b_w1); % increase in y-direction of span associated to fuselage width (tip)
y = [b_w1_fus_in b_w1_fus_out]'/2;
c_y = [cR_w cT_w1]';
le_y_1 = (x_w1R_LE - x_w1R_LE);
le_y_2 = (x_w1T_LE - x_w1R_LE);
le_y = [le_y_1 le_y_2]';
diedro = [dihedral_w dihedral_w];
[x_mesh,y_mesh,z_mesh,CA,MISC] = get_DATA_aero(NACA_foil,S,b,y,c_y,le_y,diedro,VTP,t_c);
