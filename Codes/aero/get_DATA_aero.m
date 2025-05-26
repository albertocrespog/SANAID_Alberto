function [x_mesh, y_mesh, z_mesh, CA, MISC] = get_DATA_aero(NACA_foil,S,b,y,c_y,le_y,diedro,VTP,t_c)

if NACA_foil == 1
    t_c = 0.12;
    [x,z]   = getNACA(t_c);
else
    load Bell_A821201_v4.dat
    x = Bell_A821201_v4(:,1)';
    z = Bell_A821201_v4(:,2)';
end

% Constants
f2m = 0.3048;
D2R = pi/180;
R2D = 180/pi;
                    
% MAC = S/b;
[cr, ct, xca, yca, zca, LAMc4, LAMc2, LAM, MISC] = aeroSurf_geom(S, b, y, c_y, le_y,diedro);

GEO.LAM = LAM;
GEO.LAMc2 = LAMc2;
GEO.LAMc4 = LAMc4;
ele1 = y/b*2;
ele2 = le_y/cr;
ele3 = c_y/cr;

chordDistr =    [y/b*2,le_y/cr,c_y/cr];
TR      = ct/cr;
zca     = yca*sin(diedro);
le_y    = chordDistr(:,2)*cr;    % Coordenadas x del borde de ataque de ala
c_y     = chordDistr(:,3)*cr;    % Cuerdas de las secciones del ala.
y       = chordDistr(:,1)*b/2;   % Coordenadas y de las secciones del ala
                       
n       = length(x);

% Stores aerodynamic center calculated
CA.xca = xca;
CA.yca = yca;
CA.zca = zca;

if VTP == 1
    x_mesh  = c_y*x + le_y*ones(1,n);
    y_mesh  = y*ones(1,n);
    z_mesh = c_y*z + y_mesh*tan(diedro(1));

    x_mesh  = [x_mesh];
    y_mesh  = [y_mesh];
    z_mesh  = [z_mesh];
else
    x_mesh  = c_y*x + le_y*ones(1,n);
    y_mesh  = y*ones(1,n);
    z_mesh  = c_y*z + y_mesh*tan(diedro(1));

    x_mesh  = [x_mesh(end:(-1):1,:); x_mesh];
    y_mesh  = [-y_mesh(end:(-1):1,:); y_mesh];
    z_mesh  = [z_mesh(end:(-1):1,:); z_mesh];
end
