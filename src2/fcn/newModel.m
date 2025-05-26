function model = newModel(mainDir)
    model.name              = [];
    model.conf              = 'convencional'; % convencional || canard || convencional_canard || flyWing
    model.confVert          = 'convencional'; % convencional || twin_vertical || no_vert

%% FIELD: general
    model.general.mtow      = [];
    model.general.w_w0      = [];
    model.general.W         = [];
    model.general.Sref      = [];
    model.general.h         = [];
    model.general.rhoinf    = [];
    model.general.Tinf      = [];
    model.general.pinf      = [];
    model.general.qinf      = [];
    model.general.Vinf      = [];
    model.general.Minf      = [];
    model.general.Xcg       = [];
    model.general.Vstall    = [];
    model.general.polar     = []; % [C_D0, k1, k2] C_D = C_D0 + k1*C_L + k2*C_L^2
    model.general.Ixx       = [];
    model.general.Iyy       = [];
    model.general.Izz       = [];
    model.general.Ixz       = [];
    model.general.L         = [];
    model.general.Xna       = [];
    
    
%% FIELD: propulsion
    model.propulsion.Pmax       = [];
    model.propulsion.deltaT     = [];
    model.propulsion.F_OEI      = [];
    model.propulsion.X          = [];
    model.propulsion.Y          = [];
    model.propulsion.Z          = [];
    model.propulsion.Acoef      = [];
    model.propulsion.Bcoef      = [];
    model.propulsion.Ccoef      = [];
    model.propulsion.n          = [];
    model.propulsion.Dprop      = [];
    model.propulsion.rendProp   = [];
    model.propulsion.w_R_30     = 0.1050;
    model.propulsion.w_R_60     = 0.1460;
    model.propulsion.w_R_90     = 0.0900;
    model.propulsion.nBlades    = [];
    model.propulsion.betaBlade  = 15;
    
%% FIELD: fuselaje
    
    nPoint  = 47;
    nSec    = 9;
    skip    = 2;
    fusPath             = [mainDir,'\auxFiles\fuselageExample_piaggio.txt'];
    
    model.fuselaje      = getAdimFus(fusPath, nPoint, nSec, skip);
    model.fuselaje.CLa  = [];
    
    
%     model.fuselaje.D        = fus.D;   % Fuselage maximun depth
%     model.fuselaje.W        = fus.W;   % Fuselage maximun width
%     model.fuselaje.l        = fus.l;   % Fuselage lenth
%     model.fuselaje.Sside    = fus.Sside;   % Fuselage total lateral surface
%     model.fuselaje.Sfront   = fus.Sfront;   % Fuselage front section area
%     model.fuselaje.Stop     = fus.Stop;   % Fuselage top section area
%     model.fuselaje.vol      = fus.vol;   % Fuselage volume
%     model.fuselaje.x        = fus.x;   % x Coordinates of depth and width data 
%     model.fuselaje.D_x      = fus.D_x;   % Depth evolution with x
%     model.fuselaje.W_x      = fus.W_x;   % Width evolution with x
%     model.fuselaje.S_x      = fus.S_x;   % Section area evolution with x
%     model.fuselaje.alpha_x  = fus.alpha_x;
%     model.fuselaje.CLa      = [];   % Fuselage lift slope
%     model.fuselaje.dSdX     = fus.dSdX;
%     model.fuselaje.meshData = fus.plotData;
    
    
%     model.fuselaje.Xcg      = [];

%% FIELD: ala
    model.ala.S             = [];
    model.ala.Se            = [];
    model.ala.b             = []; 
    model.ala.AR            = [];
    model.ala.TR            = []; 
    model.ala.cr            = [];
    model.ala.ct            = [];
    model.ala.t_c           = [];   % Relación entre el espesor y la cuerda del perfil
    model.ala.MAC           = [];
    model.ala.LAM           = []; 
    model.ala.LAMc2         = []; 
    model.ala.LAMc4         = [];
    model.ala.diedro        = [];
    model.ala.y             = [];
    model.ala.le_y          = [];
    model.ala.c_y           = [];
    model.ala.Xca           = [];
    model.ala.Zca           = [];
    model.ala.xca           = [];
    model.ala.yca           = [];
    model.ala.zca           = [];
    model.ala.i             = [];
    
    model.ala.CLa           = [];
    model.ala.CL0           = [];
    model.ala.CM0           = [];
    model.ala.Cla           = [];
    model.ala.eta           = [];
    
    model.ala.y0_b2         = [];
    model.ala.y1_b2         = [];
    model.ala.cm_c          = []; % Relación entre la cuerda de la superficie movil y la cuerda del ala
    
%% FIELD: vertical    
    model.vertical.S        = []; 
    model.vertical.b        = []; 
    model.vertical.AR       = [];
    model.vertical.TR       = []; 
    model.vertical.cr       = [];
    model.vertical.ct       = [];
    model.vertical.MAC      = [];
    model.vertical.LAM      = []; 
    model.vertical.LAMc2    = []; 
    model.vertical.LAMc4    = [];
    model.vertical.z        = [];
    model.vertical.le_z     = [];
    model.vertical.c_z      = [];
    model.vertical.Xca      = [];
    model.vertical.Yca      = 0;
    model.vertical.Zca      = [];
    model.vertical.xca      = [];
    model.vertical.zca      = [];
    
    model.vertical.CLa      = [];
    model.vertical.CL0      = [];
    model.vertical.CM0      = [];
    model.vertical.Cla      = [];
    model.vertical.eta      = [];
    
    model.vertical.y0_b2    = [];
    model.vertical.y1_b2    = [];
    model.vertical.cm_c     = []; % Relación entre la cuerda movil y la total 
    model.vertical.t_c      = [];

%% FIELD: horizontal
    model.horizontal.S      = []; 
    model.horizontal.Se     = [];
    model.horizontal.b      = []; 
    model.horizontal.AR     = [];
    model.horizontal.TR     = []; 
    model.horizontal.cr     = [];
    model.horizontal.ct     = [];
    model.horizontal.t_c    = [];   % Relación entre el espesor y la cuerda del perfil
    model.horizontal.MAC    = [];
    model.horizontal.LAM    = []; 
    model.horizontal.LAMc2  = []; 
    model.horizontal.LAMc4  = [];
    model.horizontal.diedro = [];
    model.horizontal.y      = [];
    model.horizontal.le_y   = [];
    model.horizontal.c_y    = [];
    model.horizontal.Xca    = [];
    model.horizontal.Zca    = [];
    model.horizontal.xca    = [];
    model.horizontal.yca    = [];
    model.horizontal.zca    = [];
    model.horizontal.i      = [];
    
    model.horizontal.CLa    = [];
    model.horizontal.CL0    = [];
    model.horizontal.CM0    = [];
    model.horizontal.Cla    = [];
    model.horizontal.eta    = [];
    
    model.horizontal.y0_b2  = [];
    model.horizontal.y1_b2  = [];
    model.horizontal.cm_c   = []; % Relación entre la cuerda de la superficie movil y la cuerda del ala
    
%% FIELD: canard
    model.canard.S          = []; 
    model.canard.Se         = [];
    model.canard.b          = []; 
    model.canard.AR         = [];
    model.canard.TR         = []; 
    model.canard.cr         = [];
    model.canard.ct         = [];
    model.canard.t_c        = [];   % Relación entre el espesor y la cuerda del perfil
    model.canard.MAC        = [];
    model.canard.LAM        = []; 
    model.canard.LAMc2      = []; 
    model.canard.LAMc4      = [];
    model.canard.diedro     = [];
    model.canard.y          = [];
    model.canard.le_y       = [];
    model.canard.c_y        = [];
    model.canard.Xca        = [];
    model.canard.Zca        = [];
    model.canard.xca        = [];
    model.canard.yca        = [];
    model.canard.zca        = [];
    model.canard.i          = [];
    
    model.canard.CLa        = [];
    model.canard.CL0        = [];
    model.canard.CM0        = [];
    model.canard.Cla        = [];
    model.canard.eta        = [];
    
    model.canard.y0_b2      = [];
    model.canard.y1_b2      = [];
    model.canard.cm_c       = []; % Relación entre la cuerda de la superficie movil y la cuerda del ala
    
    
%% FIELD: derivadas
    % Generales
    model.derivadas.CL              = [];
    model.derivadas.CM              = [];
    model.derivadas.CD              = [];
    % Longitudinales
    model.derivadas.CL_a            = [];
    model.derivadas.CD_a            = [];
    model.derivadas.CM_a            = [];
    model.derivadas.CL_q            = [];
    model.derivadas.CD_q            = [];
    model.derivadas.CM_q            = [];
    model.derivadas.CL_u            = [];
    model.derivadas.CD_u            = [];
    model.derivadas.CM_u            = [];
    model.derivadas.CL_alphaDot     = [];
    model.derivadas.CD_alphaDot     = [];
    model.derivadas.CM_alphaDot     = [];
    % Propulsivas
    model.derivadas.CM_Ta           = [];
    model.derivadas.CT_x1           = [];
    model.derivadas.CT_xu           = [];
    model.derivadas.CT_xa           = [];
    model.derivadas.CM_Tu           = [];
    model.derivadas.CM_T1           = [];
    model.derivadas.Cy_Tbeta        = [];
    model.derivadas.Cn_Tbeta        = [];
    % Laterales-direccionales
    model.derivadas.Cy_beta         = [];
    model.derivadas.Cl_beta         = [];
    model.derivadas.Cn_beta         = [];
    model.derivadas.Cy_p            = [];
    model.derivadas.Cl_p            = [];
    model.derivadas.Cn_p            = [];
    model.derivadas.Cy_r            = [];
    model.derivadas.Cl_r            = [];
    model.derivadas.Cn_r            = [];
    model.derivadas.Cy_betaDot      = [];
    model.derivadas.Cl_betaDot      = [];
    model.derivadas.Cn_betaDot      = [];
    % Control longitudinal
    model.derivadas.CD_de           = [];
    model.derivadas.CL_de           = [];
    model.derivadas.CM_de           = [];
    model.derivadas.CD_dc           = [];
    model.derivadas.CL_dc           = [];
    model.derivadas.CM_dc           = [];
    % Control lateral
    model.derivadas.Cy_da           = [];
    model.derivadas.Cl_da           = [];
    model.derivadas.Cn_da           = [];
    model.derivadas.Cy_dr           = [];
    model.derivadas.Cl_dr           = [];
    model.derivadas.Cn_dr           = [];
    
    