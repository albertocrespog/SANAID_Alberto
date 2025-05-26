%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   19 de Julio de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
% Proporciona las incidencias de 2 superficies fijadas las posiciones de
% estas y la incicidencia de una de ellas en caso de contar con 3
% superficies horizontales

function [i_out] = pre_equilibrio(model)

%% ADQUISICIÓN DE DATOS DEL MODELO
% Datos generales
S_ref   = model.general.Sref;
W       = model.genera.mtow*9.8065;
q       = model.general.q;
X_cg    = model.general.Xcg;

% Datos de ala
X_w     = model.ala.Xca;
S_w     = model.ala.S;
CLa_w   = model.ala.CLa;
CL0_w   = model.ala.CL0;
CM_w    = model.ala.CM;
c_w     = model.ala.c;

%% RESOLUCION
% La solución es obtenida mediante resolución simbólica del sistema de 
% ecuaciones indicadas a continuación que no son más que las ecuaciones de
% equilibrio de fuerzas y momentos.

%     syms i_t_var i_c_var i_w_var
%     
%     eq_F =  CL0_w + eta_c*S_c/S_ref*CL0_c + eta_t*S_t/S_ref*CL0_t + ...
%             CLa_w*i_w_var + eta_c*S_c/S_ref*CLa_c*(i_c_var + E0_c) + ...
%             eta_t*S_t/S_ref*CLa_t*(i_t_var - E0_t) == W/q/S_ref;
%     
%     eq_M =  eta_c*S_c/S_ref*c_c/c_w*CM_c + eta_c*S_c/S_ref/c_w*(X_cg-X_c)...
%             *(CL0_c + CLa_c*(i_c_var + E0_c)) + CM_w + S_w/S_ref/c_w*...
%             (X_cg-X_w)*(CL0_w+CLa_w*i_w_var)+eta_t*S_t/S_ref*c_t/c_w*CM_t...
%             + eta_t*S_t/S_ref*(X_cg-X_t)*(CL0_t+CLa_t*(i_t_var - E0_t))==0;

% Este método tiene el inconveniente de ser un poco lento, sobre todo la
% primera vez que se ejecuta ya que no esta iniciado el modulo Symbolic
% Math Toolbox

switch model.conf  % Switch entre las posibles configuraciones
    case 'convencional' % Ala + tail
        X_t     = model.horizontal.Xca;
        S_t     = model.horizontal.S;
        CLa_t   = model.horizontal.CLa;
        CL0_t   = model.horizontal.CL0;
        CM_t    = model.horizontal.CM;
        c_t     = model.horizontal.c;
        eta_t   = model.horizontal.eta;
        E0_t    = 0;
        
        syms i_w i_t;
        eq_F    = S_w/S_ref*CL0_w + S_w/S_ref*CLa_w*i_w + (CL0_t*S_t*eta_t)/S_ref + ...
                (CLa_t*S_t*eta_t*(i_t - E0_t))/S_ref == W/(S_ref*q);
        
        eq_M    = S_w/S_ref*CM_w + (S_w*(X_cg - X_w)*(CL0_w + CLa_w*i_w))/(S_ref*c_w)...
                + (S_t*eta_t*(X_cg - X_t)*(CL0_t + CLa_t*(i_t - E0_t)))/...
                S_ref + (CM_t*S_t*c_t*eta_t)/(S_ref*c_w) == 0;
        
        sol = solve(eq_F,eq_M);
        i_out(1) = double(sol.i_w);
        i_out(2) = double(sol.i_t);
        
    case 'canard' % Ala + canard
        
        X_c     = model.canard.Xca;
        S_c     = model.canard.S;
        CLa_c   = model.canard.CLa;
        CL0_c   = model.canard.CL0;
        CM_c    = model.canard.CM;
        c_c     = model.canard.c;
        eta_c   = model.canard.eta;
        E0_c    = 0; 
        
        syms i_w i_c;
        eq_F    = S_w/S_ref*CL0_w + S_w/S_ref*CLa_w*i_w + (CL0_c*S_c*eta_c)/S_ref + ...
                (CLa_c*S_c*eta_c*(E0_c + i_c))/S_ref == W/(S_ref*q);
        
        eq_M    = S_w/S_ref*CM_w + (X_cg - X_w)*(CL0_w + CLa_w*i_w)*S_w/(S_ref*c_w)...
                + (CM_c*S_c*c_c*eta_c)/(S_ref*c_w) + (S_c*eta_c*(CL0_c +... 
                CLa_c*(E0_c + i_c))*(X_cg - X_c))/(S_ref*c_w) == 0;
        
        sol = solve(eq_F,eq_M);
        i_out(1) = double(sol.i_w);
        i_out(2) = double(sol.i_c);
        
    case 'convencional_canard' % Ala + tail + canard
        
        X_t     = model.horizontal.Xca;
        S_t     = model.horizontal.S;
        CLa_t   = model.horizontal.CLa;
        CL0_t   = model.horizontal.CL0;
        CM_t    = model.horizontal.CM;
        c_t     = model.horizontal.c;
        eta_t   = model.horizontal.eta;
        E0_t    = 0;
        
        X_c     = model.canard.Xca;
        S_c     = model.canard.S;
        CLa_c   = model.canard.CLa;
        CL0_c   = model.canard.CL0;
        CM_c    = model.canard.CM;
        c_c     = model.canard.c;
        eta_c   = model.canard.eta;
        E0_c    = 0;
        
        syms i_t i_c i_w
        eq_F =  S_w/S_ref*CL0_w + eta_c*S_c/S_ref*CL0_c + eta_t*S_t/S_ref*CL0_t + ...
        S_w/S_ref*CLa_w*i_w + eta_c*S_c/S_ref*CLa_c*(i_c + E0_c) + ...
        eta_t*S_t/S_ref*CLa_t*(i_t - E0_t) == W/q/S_ref;
        
        eq_M =  eta_c*S_c/S_ref*c_c/c_w*CM_c + eta_c*S_c/S_ref/c_w*(X_cg-X_c)*...
        (CL0_c + CLa_c*(i_c + E0_c)) + S_w/S_ref*CM_w + S_w/S_ref/c_w*...
        (X_cg-X_w)*(CL0_w+CLa_w*i_w)+eta_t*S_t/S_ref*c_t/c_w*CM_t...
        + eta_t*S_t/S_ref*(X_cg-X_t)*(CL0_t+CLa_t*(i_t - E0_t))==0;
  
        switch model.nfo.i_fix     % Switch entre la incidencia fijada
            case 'w'    % Esta fijada la incidencia del ala
                i_w_val = model.ala.i;   
                
                eq_F = subs(eq_F,i_w,i_w_val);
                eq_M = subs(eq_M,i_w,i_w_val);
        
                sol = solve(eq_F,eq_M);
                i_out(1) = double(sol.i_t);
                i_out(2) = double(sol.i_c);
            
            case 'c'    % Esta fijada la incidencia del canard    
                i_c_val = model.canard.i;           
                
                eq_F = subs(eq_F,i_c,i_c_val);
                eq_M = subs(eq_M,i_c,i_c_val);
            
                sol = solve(eq_F,eq_M);
                i_out(1) = double(sol.i_w);
                i_out(2) = double(sol.i_t);
            
            case 't'    % Esta fijada la incidencia de la cola
                i_t_val = model.horizontal.i;   
                
                eq_F = subs(eq_F,i_t,i_t_val);
                eq_M = subs(eq_M,i_t,i_t_val);
                       
                sol = solve(eq_F,eq_M);
                i_out(1) = double(sol.i_w);
                i_out(2) = double(sol.i_c);

        end
end
end

% syms i_t_var i_c_var i_w_var S_ref S_w S_t S_c X_w X_t X_c X_cg c_w CLa_w CLa_t CLa_c eta_t eta_c W q CL0_w CL0_c CL0_t CM_w CM_c CM_t E0_c E0_t c_c c_t