function result = trim_longitudinal(modelo, Xcg_evol, W_W0_evol, resolMode, fix_surf)

n       = length(Xcg_evol);


%% Datos generales
g       = 9.8065;
e       = 0.95; % Oswald Efficiency
W0      = modelo.general.mtow*g;
Sref    = modelo.general.Sref;
q       = modelo.general.qinf;
V       = modelo.general.Vinf;


%% Datos ala

CLa_a   = modelo.ala.CLa;
rend_a  = modelo.ala.eta;
CL0_a   = modelo.ala.CL0;
CM_a    = modelo.ala.CM0;
AR_a    = modelo.ala.AR;
Sa      = modelo.ala.S;
ba      = modelo.ala.b;
TR_a    = modelo.ala.TR;
Xa      = modelo.ala.Xca;
i_a     = modelo.ala.i;
c       = modelo.ala.cr;
LAMc4a  = modelo.ala.LAMc4;
za      = modelo.ala.Zca;

K       = 1/(pi*AR_a*e); 


%% Datos superficies estabilizadoras

switch modelo.conf
    case 'canard'
        CLa_h               = 0;
        CL0_h               = 0;
        CM_h                = 0;
        rend_h              = 0;
        Sh                  = 0;
        i_h                 = 0;
        Xh                  = 0;
        downw               = 0;
        alphaCL_alphaCl_h   = 0;
        Clde_Cldetheory_h   = 0;
        Cldetheory_h        = 0;
        Kb_h                = 0;
        Cla_h               = 0;
        
        CLa_c   = modelo.canard.CLa;
        rend_c  = modelo.canard.eta;
        CL0_c   = modelo.canard.CL0;
        CM_c    = modelo.canard.CM0;
        Sc      = modelo.canard.S;
        Xc      = modelo.canard.Xca;
        i_c     = modelo.canard.i;
        TR_c    = modelo.canard.TR;
        Cla_c   = modelo.canard.Cla;
        y1_b2_c = 2*modelo.canard.y0_b2;
        y2_b2_c = 2*modelo.canard.y1_b2;
        cf_c_c  = modelo.canard.cm_c;
        t_c_c   = modelo.canard.t_c;
        AR_c    = modelo.canard.AR;
        
        Kb_c                = Kb_calc(y1_b2_c, y2_b2_c, TR_c);
        Clde_Cldetheory_c   = Cldelta_Cldeltatheory_calc(cf_c_c, Cla_c); 
        Cldetheory_c        = Cldeltatheory_calc(cf_c_c, t_c_c);
        alphaCL_alphaCl_c   = alphaCL_alphaCl_calc(cf_c_c, AR_c);
        upw                 = upwash_calc(AR_a, Xa, Xc, c);
        
        Kh  = 0;
        Kc = 1/(pi*AR_c*e);
        
    case 'convencional'
        CLa_c               = 0;
        CL0_c               = 0;
        CM_c                = 0;
        rend_c              = 0;
        Sc                  = 0;
        i_c                 = 0;
        Xc                  = 0;
        upw                 = 0;
        alphaCL_alphaCl_c   = 0;
        Clde_Cldetheory_c   = 0;
        Cldetheory_c        = 0;
        Kb_c                = 0;
        Cla_c               = 0;
        
        CLa_h   = modelo.horizontal.CLa;
        rend_h  = modelo.horizontal.eta;
        CL0_h   = modelo.horizontal.CL0;
        CM_h    = modelo.horizontal.CM0;
        Sh      = modelo.horizontal.S;
        Xh      = modelo.horizontal.Xca;
        i_h     = modelo.horizontal.i;
        TR_h    = modelo.horizontal.TR;
        Cla_h   = modelo.horizontal.Cla;
        y1_b2_h = modelo.horizontal.y0_b2;
        y2_b2_h = modelo.horizontal.y1_b2;
        cf_c_h  = modelo.horizontal.cm_c;
        t_c_h   = modelo.horizontal.t_c;
        AR_h    = modelo.horizontal.AR;
        zh      = modelo.horizontal.Zca;
        
        Kb_h                = Kb_calc(y1_b2_h, y2_b2_h, TR_h);
        Clde_Cldetheory_h   = Cldelta_Cldeltatheory_calc(cf_c_h, Cla_h); 
        Cldetheory_h        = Cldeltatheory_calc(cf_c_h, t_c_h);
        alphaCL_alphaCl_h   = alphaCL_alphaCl_calc(cf_c_h ,AR_h);
        downw               = downwash_calc(Xa, Xh, AR_a, TR_a, LAMc4a, ba, zh - za);
        
        Kc  = 0;
        Kh   = 1/(pi*AR_h*e);
        
        
    case 'flyWing'
        CLa_c               = 0;
        CL0_c               = 0;
        CM_c                = 0;
        rend_c              = 0;
        Sc                  = 0;
        i_c                 = 0;
        Xc                  = 0;
        upw                 = 0;
        alphaCL_alphaCl_c   = 0;
        Clde_Cldetheory_c   = 0;
        Cldetheory_c        = 0;
        Kb_c                = 0;
        Cla_c               = 0;
        
        % Virtual Horizontal
        CLa_h   = modelo.horizontal.CLa;
        rend_h  = modelo.horizontal.eta;
        CL0_h   = modelo.horizontal.CL0;
        CM_h    = modelo.horizontal.CM0;
        Sh      = modelo.horizontal.S;
        Xh      = modelo.horizontal.Xca;
        i_h     = modelo.horizontal.i;
        TR_h    = modelo.horizontal.TR;
        Cla_h   = modelo.horizontal.Cla;
        y1_b2_h = modelo.horizontal.y0_b2;
        y2_b2_h = modelo.horizontal.y1_b2;
        cf_c_h  = modelo.horizontal.cm_c;
        t_c_h   = modelo.horizontal.t_c;
        AR_h    = modelo.horizontal.AR;
        zh      = modelo.horizontal.Zca;
        
        Kb_h                = Kb_calc(y1_b2_h, y2_b2_h, TR_h);
        Clde_Cldetheory_h   = Cldelta_Cldeltatheory_calc(cf_c_h, Cla_h); 
        Cldetheory_h        = Cldeltatheory_calc(cf_c_h, t_c_h);
        alphaCL_alphaCl_h   = alphaCL_alphaCl_calc(cf_c_h ,AR_h);
        downw               = downwash_calc(Xa, Xh, AR_a, TR_a, LAMc4a, ba, zh - za);
        
        Kc  = 0;
        Kh   = 1/(pi*AR_h*e);
        
    case 'convencional_canard'
        CLa_h   = modelo.horizontal.CLa;
        rend_h  = modelo.horizontal.eta;
        CL0_h   = modelo.horizontal.CL0;
        CM_h    = modelo.horizontal.CM0;
        Sh      = modelo.horizontal.S;
        Xh      = modelo.horizontal.Xca;
        i_h     = modelo.horizontal.i;
        TR_h    = modelo.horizontal.TR;
        Cla_h   = modelo.horizontal.Cla;
        y1_b2_h = modelo.horizontal.y0_b2;
        y2_b2_h = modelo.horizontal.y1_b2;
        cf_c_h  = modelo.horizontal.cm_c;
        t_c_h   = modelo.horizontal.t_c;
        AR_h    = modelo.horizontal.AR;
        zh      = modelo.horizontal.Zca;
        
        Kb_h                = Kb_calc(y1_b2_h, y2_b2_h, TR_h);
        Clde_Cldetheory_h   = Cldelta_Cldeltatheory_calc(cf_c_h, Cla_h); 
        Cldetheory_h        = Cldeltatheory_calc(cf_c_h, t_c_h);
        alphaCL_alphaCl_h   = alphaCL_alphaCl_calc(cf_c_h ,AR_h);
        downw               = downwash_calc(Xa, Xh, AR_a, TR_a, LAMc4a, ba, zh - za);
        
        CLa_c   = modelo.canard.CLa;
        rend_c  = modelo.canard.eta;
        CL0_c   = modelo.canard.CL0;
        CM_c    = modelo.canard.CM0;
        Sc      = modelo.canard.S;
        Xc      = modelo.canard.Xca;
        i_c     = modelo.canard.i;
        TR_c    = modelo.canard.TR;
        Cla_c   = modelo.canard.Cla;
        y1_b2_c = 2*modelo.canard.y0_b2;
        y2_b2_c = 2*modelo.canard.y1_b2;
        cf_c_c  = modelo.canard.cm_c;
        t_c_c   = modelo.canard.t_c;
        AR_c    = modelo.canard.AR;
        
        Kb_c                = Kb_calc(y1_b2_c, y2_b2_c, TR_c);
        Clde_Cldetheory_c   = Cldelta_Cldeltatheory_calc(cf_c_c, Cla_c); 
        Cldetheory_c        = Cldeltatheory_calc(cf_c_c, t_c_c);
        alphaCL_alphaCl_c   = alphaCL_alphaCl_calc(cf_c_c, AR_c);
        upw                 = upwash_calc(AR_a, Xa, Xc, c);

        switch fix_surf
            case 'canard'
                Kc  = 0;
                Kh   = 1/(pi*AR_h*e);
            case 'horizontal'
                Kh  = 0;
                Kc = 1/(pi*AR_c*e);
        end
end



%{        
if isempty(modelo.horizontal.S)
    CLa_h               = 0;
    CL0_h               = 0;
    CM_h                = 0;
    rend_h              = 0;
    Sh                  = 0;
    i_h                 = 0;
    Xh                  = 0;
    downw               = 0;
    alphaCL_alphaCl_h   = 0;
    Clde_Cldetheory_h   = 0;
    Cldetheory_h        = 0;
    Kb_h                = 0;
    Cla_h               = 0;
else
    CLa_h   = modelo.horizontal.CLa;
    rend_h  = modelo.horizontal.eta;
    CL0_h   = modelo.horizontal.CL0;
    CM_h    = modelo.horizontal.CM0;
    Sh      = modelo.horizontal.S;
    Xh      = modelo.horizontal.Xca;
    i_h     = modelo.horizontal.i;
    TR_h    = modelo.horizontal.TR;
    Cla_h   = modelo.horizontal.Cla;
    y1_b2_h = modelo.horizontal.y0_b2;
    y2_b2_h = modelo.horizontal.y1_b2;
    cf_c_h  = modelo.horizontal.cm_c;
    t_c_h   = modelo.horizontal.t_c;
    AR_h    = modelo.horizontal.AR;
    zh      = modelo.horizontal.Zca;
        
    Kb_h                = Kb_calc(y1_b2_h, y2_b2_h, TR_h);
    Clde_Cldetheory_h   = Cldelta_Cldeltatheory_calc(cf_c_h, Cla_h); 
    Cldetheory_h        = Cldeltatheory_calc(cf_c_h, t_c_h);
    alphaCL_alphaCl_h   = alphaCL_alphaCl_calc(cf_c_h ,AR_h);
    downw               = downwash_calc(Xa, Xh, AR_a, TR_a, LAMc4a, ba, zh - za);
end

if isempty(modelo.canard.S)   
    CLa_c               = 0;
    CL0_c               = 0;
    CM_c                = 0;
    rend_c              = 0;
    Sc                  = 0;
    i_c                 = 0;
    Xc                  = 0;
    upw                 = 0;
    alphaCL_alphaCl_c   = 0;
    Clde_Cldetheory_c   = 0;
    Cldetheory_c        = 0;
    Kb_c                = 0;
    Cla_c               = 0;
else
    CLa_c   = modelo.canard.CLa;
    rend_c  = modelo.canard.eta;
    CL0_c   = modelo.canard.CL0;
    CM_c    = modelo.canard.CM0;
    Sc      = modelo.canard.S;
    Xc      = modelo.canard.Xca;
    i_c     = modelo.canard.i;
    TR_c    = modelo.canard.TR;
    Cla_c   = modelo.canard.Cla;
    y1_b2_c = 2*modelo.canard.y0_b2;
    y2_b2_c = 2*modelo.canard.y1_b2;
    cf_c_c  = modelo.canard.cm_c;
    t_c_c   = modelo.canard.t_c;
    AR_c    = modelo.canard.AR;
        
    Kb_c                = Kb_calc(y1_b2_c, y2_b2_c, TR_c);
    Clde_Cldetheory_c   = Cldelta_Cldeltatheory_calc(cf_c_c, Cla_c); 
    Cldetheory_c        = Cldeltatheory_calc(cf_c_c, t_c_c);
    alphaCL_alphaCl_c   = alphaCL_alphaCl_calc(cf_c_c, AR_c);
    upw                 = upwash_calc(AR_a, Xa, Xc, c);
end
%}

P_total = modelo.propulsion.Pmax;
z_motor = modelo.propulsion.Z;

l_fus   = modelo.general.L;



%% Ala
% CLi_a   = rend_a*Sa/Sref*CLa_a;
% CMi_a   = - rend_a*Sa/Sref*CLa_a*(Xa-Xcg)/c;

%% Canard
% CLi_c       = rend_c*Sc/Sref*CLa_c;
% CMi_c       = - rend_c*Sc/Sref*CLa_c*(Xc-Xcg)/c;
% alpha_de_c  = Kb_c*Clde_Cldetheory_c*Cldetheory_c/Cla_c*alphaCL_alphaCl_c ;
% CLde_c      = alpha_de_c*CLi_c;
% CMde_c      = alpha_de_c*CMi_c;


%% Horizontal
% CLi_h       = rend_h*Sh/Sref*CLa_h;
% CMi_h       = - rend_h*Sh/Sref*CLa_h*(Xh-Xcg)/c;
% alpha_de_h  = Kb_h*Clde_Cldetheory_h*Cldetheory_h/Cla_h*alphaCL_alphaCl_h;
% CLde_h      = alpha_de_h*CLi_h;
% CMde_h      = alpha_de_h*CMi_h;

%% Fuselaje

if strcmp(resolMode,'basic')
    CLi_f = 0;
    CMi_f = 0;
    P_total = 0;
    z_motor = 0;
else
    CLi_f = Sa/Sref*CLa_f;
    CMi_f = Sa/Sref*CMa_f;    
end



%% Caracteristicas del avion

% CLa = CLi_a + CLi_c*(1 + upw) + CLi_h*(1-downw) + CLi_f;
% 
% CMa = CMi_a + CMi_h + CMi_c + CMi_f;
% 
% CL0 = CLi_a*i_a + rend_a*Sa/Sref*CL0_a + CLi_c*(i_c + e0c) + rend_c*Sc/Sref*CL0_c + CLi_h*(i_h - e0h) + rend_h*Sh/Sref*CL0_h;
% 
% CM0 =   Sa/Sref*rend_a*CM_a + (CLi_a*i_a + rend_a*Sa/Sref*CL0_a)*(Xw - Xcg)/c + Sc/Sref*rend_c*CM_c + (CLi_c*(i_c + upw) + rend_c*Sc/Sref*CL0_c)*(Xc - Xcg)/c...
%         + Sh/Sref*rend_h*CM_h + (CLi_h*(i_h - downw) + rend_h*Sh/Sref*CL0_h)*(Xh - Xcg)/c;  


% 2*MTOW*W_W0/rho/Sref/V^2 = CL0 + CLa*alpha + CLdc*de_c + CLde*de_h + P_total/V*sen(alpha_motor);
% 0 = CMa*alpha + CMde_c*de_c + CMde_h*de_h + (CL0 + CLa*alpha)*(Xna - Xcg)/c + CLde_c*de_c*(Xc - Xcg)/c + CLde_h*de_h*(Xh - Xcg)/c - P_total/V*z_motor

result = zeros(7,n);
CLi_a   = rend_a*Sa/Sref*CLa_a;

e0h  = 2*CL0_a/(pi*AR_a);
e0c  = e0h;
if Sh == 0
    CLde_h      = 0;
    alpha_de_h  = 0;
    CLi_h       = 0;
else
    CLi_h       = rend_h*Sh/Sref*CLa_h;
    alpha_de_h  = Kb_h*Clde_Cldetheory_h*Cldetheory_h/Cla_h*alphaCL_alphaCl_h;
    CLde_h      = alpha_de_h*CLi_h;
end
if Sc == 0
    CLde_c      = 0;
    alpha_de_c  = 0;
    CLi_c       = 0;
else
    CLi_c       = rend_c*Sc/Sref*CLa_c;
    alpha_de_c  = Kb_c*Clde_Cldetheory_c*Cldetheory_c/Cla_c*alphaCL_alphaCl_c;
    CLde_c      = alpha_de_c*CLi_c;
end

if strcmp(modelo.conf,'flyWing')
    CLa = CLi_a + CLi_f;
    CL0 = CLi_a*i_a + rend_a*Sa/Sref*CL0_a;
else
    CLa = CLi_a + CLi_c*(1 + upw) + CLi_h*(1-downw) + CLi_f;
    CL0 = CLi_a*i_a + rend_a*Sa/Sref*CL0_a + CLi_c*(i_c + e0c) + rend_c*Sc/Sref*CL0_c + CLi_h*(i_h - e0h) + rend_h*Sh/Sref*CL0_h;
end
result(4,:) = CL0;
result(6,:) = CLa;

for kk = 1:n
    Xcg     = l_fus*Xcg_evol(kk);
    W       = W0*W_W0_evol(kk);
    % Coeficiente de momentos ala
    CMi_a   = - rend_a*Sa/Sref*CLa_a*(Xa-Xcg)/c;
    
    % Coeficiente de momentos del horizontal
    CMi_h       = - rend_h*Sh/Sref*CLa_h*(Xh-Xcg)/c;
    CMde_h      = alpha_de_h*CMi_h;
    
    % Coeficiente de momentos del canard
    CMi_c       = - rend_c*Sc/Sref*CLa_c*(Xc-Xcg)/c;
    CMde_c      = alpha_de_c*CMi_c;
    
    % Coeficiente de momentos de la aeronave
    if strcmp(modelo.conf,'flyWing')
        CMa = CMi_a;
    else
        CMa = CMi_a + CMi_h + CMi_c + CMi_f;
    end
    result(8,kk) = CMa;
    if isempty(i_c)
        i_c = 0;
    end
    if isempty(i_h)
        i_h = 0;
    end
    if isempty(i_a)
        i_a = 0;
    end
    if strcmp(modelo.conf,'flyWing')
        CM0 = Sa/Sref*rend_a*CM_a + (CLi_a*i_a + rend_a*Sa/Sref*CL0_a)*(Xcg - Xa)/c;
    else
        CM0 = Sa/Sref*rend_a*CM_a + (CLi_a*i_a + rend_a*Sa/Sref*CL0_a)*(Xcg - Xa)/c + Sc/Sref*rend_c*CM_c + (CLi_c*(i_c + e0c)...
            + rend_c*Sc/Sref*CL0_c)*(Xcg - Xc)/c + Sh/Sref*rend_h*CM_h + (CLi_h*(i_h - e0h) + rend_h*Sh/Sref*CL0_h)*(Xcg - Xh)/c;
    end
    result(5,kk) = CM0;
    switch modelo.conf
        case 'convencional'
            A = [CLa CLde_h; CMa CMde_h];
            b = [W/q/Sref - CL0; P_total/V*z_motor - CM0];
            result(7,kk) = CLde_h;
            result(9,kk) = CMde_h;
            
        case 'flyWing'
            A = [CLa CLde_h; CMa CMde_h];
            b = [W/q/Sref - CL0; P_total/V*z_motor - CM0];
            result(7,kk) = CLde_h;
            result(9,kk) = CMde_h;
            
        case 'canard'
            A = [CLa CLde_c; CMa CMde_c];
            b = [W/q/Sref - CL0; P_total/V*z_motor - CM0];
            result(7,kk) = CLde_c;
            result(9,kk) = CMde_c;
            
        case 'convencional_canard'
            switch fix_surf
                case 'canard'
                    A = [CLa CLde_h; CMa CMde_h];
                    b = [2*W/q/Sref - CL0; (P_total/V*z_motor)/(q*Sref*c) - CM0];
                    result(7,kk) = CLde_h;
                    result(9,kk) = CMde_h;
                    
                case 'horizontal'
                    A = [CLa CLde_c; CMa CMde_c];
                    b = [W/q/Sref - CL0; (P_total/V*z_motor)/(q*Sref*c) - CM0];
                    result(7,kk) = CLde_c;
                    result(9,kk) = CMde_c;
            end
    end
    result(1:2,kk) = A\b;
    
end
CDiTrim = K*(CLa*result(1,:)).^2 + rend_h*Sh/Sref*Kh*(CLde_h*result(2,:)).^2 + rend_c*Sc/Sref*Kc*(CLde_c*result(2,:)).^2;
result(3,:) = CDiTrim;

end