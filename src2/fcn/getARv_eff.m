function ARv_eff = getARv_eff(S_v, b_v, TR_v, S_h, dx_vle2hac, dz_fus2hac, Dfus, type)

AR_v = b_v^2/S_v;
MAC_v = S_v/b_v;

switch type
    case 'convencional'
        AvB_Av      = AvB_Av_calc(b_v, Dfus, TR_v);
        KH          = KH_calc(S_h, S_v);
        AvHB_AvB    = AvHB_AvB_calc(dx_vle2hac, b_v, dz_fus2hac, MAC_v);
        ARv_eff     = AvB_Av*AR_v*(1+KH*(AvHB_AvB-1));        
    case 'twin_vertical'
        % TODO: Implement real twin vertical AR effec estimation
        AvB_Av      = AvB_Av_calc(b_v, Dfus, TR_v);
        KH          = KH_calc(S_h, S_v);
        AvHB_AvB    = AvHB_AvB_calc(dx_vle2hac, b_v, dz_fus2hac, MAC_v);
        ARv_eff     = AvB_Av*AR_v*(1+KH*(AvHB_AvB-1));    
    case 'v_tail'
end

% b_v:      Envergadura del vertical
%
% Dfus:     Diametro medio del fuselaje en donde se encuentra situado el
%           vertical.
%
% TR_v:     Estrechamiento del vertical.
%
% ARv:      Relación de aspecto del vertical.
%
% S_h:       Superficie del estabilizador horizontal.
%
% Sv:       Superficie del estabilizador vertical.
%
% dx_vle2hac:    Distancia en el eje x entre el centro aerodinamico del estabilizador horizontal y el borde de ataque del estabilizador vertical.
%
% dz_fus2hac:    Distancia en el eje z entre el plano del estabilizador horizontal y la linea
%           central del fuselaje.
%
% AvB_Av:   Ratio of the vertical tail aspect ratio in the presence of the
%           fuselage to that of the isolated vertical tail.
%
% AvHB_AvB: Ratio of the vertical tail aspect ratio in the presence of the
%           fuselage and the horizontal tail to that in the presence of the fuselage alone.
%
% KH:       Factor that accounts for the relative size of the horizontal and the
%           vertical tail.
