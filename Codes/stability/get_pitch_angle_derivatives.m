function    Stab_Der = get_pitch_angle_derivatives(Stab_Der, TRIM_RESULTS)
CL = Stab_Der.CL;
% Euler angle
% theta1 = TRIM_RESULTS.trim_alpha;
theta1 = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PITCH ANGLE DERIVATES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CZ_teta = -CL*sin(theta1); % Pamadi 4.440
    CX_teta = -CL*cos(theta1); % Pamadi 4.434
    CM_teta = 0;
    
    %% Assumption!!! NEED TO CHECK!!
    CD_teta = -CX_teta;
    CL_teta = -CZ_teta;
    
    Stab_Der.CX_teta=CX_teta;
    Stab_Der.CD_teta=CD_teta;
    Stab_Der.CZ_teta=CZ_teta;
    Stab_Der.CL_teta=CL_teta;
    Stab_Der.CM_teta=CM_teta;
end