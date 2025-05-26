function OUTPUT_read_XLSX = initialilize_segments_performance_mission(SEGMENTS,OUTPUT_read_XLSX)
% Taxy
OUTPUT_read_XLSX.IPP_flags.temp_local_taxy  = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.h_inicial_taxy   = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.P_local_taxy     = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.delta_T_taxy     = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.V_taxy           = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.t_taxy           = zeros(1,length(SEGMENTS));
% Take Off
OUTPUT_read_XLSX.IPP_flags.temp_local_TO    = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.h_inicial_TO     = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.P_local_TO       = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.mu_TO            = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.h_obstacle_TO    = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.gamma_climb_TO   = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.delta_T_TO       = zeros(1,length(SEGMENTS));
% Climb
OUTPUT_read_XLSX.IPP_flags.h_inicial_cl     = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.h_final_cl       = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.gamma_cl         = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.Mach_cl          = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.TAS_cl           = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.EAS_cl           = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.delta_T_cl       = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.V_ini_cl         = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.V_fin_cl         = zeros(1,length(SEGMENTS));
% VTOL Climb
OUTPUT_read_XLSX.IPP_flags.h_inicial_vtcl   = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.h_final_vtcl     = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.t_hover          = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.delta_T_vtcl     = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.vclimb_vtcl      = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.mbat_vtcl        = zeros(1,length(SEGMENTS));
% Cruise
OUTPUT_read_XLSX.IPP_flags.h_inicial_cr     = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.dist_final_cr    = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.V_cr             = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.delta_T_cr       = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.V_ini_cr         = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.V_fin_cr         = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.fuel_cr          = zeros(1,length(SEGMENTS)); % - [kg] % 7: COMBUSTIBLE A QUEMAR
OUTPUT_read_XLSX.IPP_flags.Cd0_cr           = zeros(1,length(SEGMENTS)); % - [] % 8: CDO = F(M)  CD0: CD = CD0 + K1*CL^2 - K2*CL
OUTPUT_read_XLSX.IPP_flags.k1_cr            = zeros(1,length(SEGMENTS)); % - []% 9: K1 = F(M)  K1: CD = CD0 + K1*CL^2 - K2*CL
OUTPUT_read_XLSX.IPP_flags.k2_cr            = zeros(1,length(SEGMENTS)); % - [] % 10: K2 = F(M) K2: CD = CD0 + K1*CL^2 - K2*CL
OUTPUT_read_XLSX.IPP_flags.mbat             = zeros(1,length(SEGMENTS));
% Load Deployment
OUTPUT_read_XLSX.IPP_flags.carga_loadep     = zeros(1,length(SEGMENTS));
% Turn
OUTPUT_read_XLSX.IPP_flags.h_inicial_tr     = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.t_final_tr       = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.V_turn           = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.delta_T_tr       = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.phi_tr           = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.V_psi            = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.n_tr             = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.R_tr             = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.h_inicial_tr_wt  = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.t_final_tr_wt    = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.V_turn_wt        = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.delta_T_tr_wt    = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.phi_tr_wt        = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.V_psi_wt         = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.n_tr_wt          = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.R_tr_wt          = zeros(1,length(SEGMENTS));
% Descent
OUTPUT_read_XLSX.IPP_flags.h_inicial_d      = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.h_final_d        = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.gamma_d          = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.Mach_d           = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.EAS_d            = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.TAS_d            = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.delta_T_d        = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.V_ini_d          = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.V_fin_d          = zeros(1,length(SEGMENTS));
% VTOL Descent
OUTPUT_read_XLSX.IPP_flags.h_inicial_vtd    = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.h_final_vtd      = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.vdes_vtd         = zeros(1,length(SEGMENTS));
% Landing
OUTPUT_read_XLSX.IPP_flags.temp_local_LND   = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.h_inicial_LND    = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.P_local_LND      = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.mu_LND           = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.delta_T_LND      = zeros(1,length(SEGMENTS));
OUTPUT_read_XLSX.IPP_flags.t_brake          = zeros(1,length(SEGMENTS));
% Dummy
OUTPUT_read_XLSX.IPP_flags.dummy            = zeros(1,3);