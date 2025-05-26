function [case_AC OUTPUT_read_XLSX] = Initial_prompt_advanced

%% Selects Automatization process

% Prompt aircraft Selection
% Defines de aircraft to be analyzed
% case_AC = 1 - EMERGENTIA 1:1 case_AC = 2 - EMERGENTIA 1:2 case_AC = 3 - PEPIÑOXXL case_AC = 4 - Comercial AC case_AC = 5 - WIG case_AC = 6 - CERVERA case_AC = 7 - Existing AC - B747       % case_AC = 8 - Tamiz
prompt = sprintf('ENTER AICRAFT SELECTION:\n1 - EMERGENTIA 100%% SCALE\n2 - EMERGENTIA MANUFACTURED%% SCALE\n3 - PEPIÑOXXL\n4 - Comertial AC\n5 - WIG\n6 - CERVERA\n7 - Existing AC,\n8 - Tamiz\n9 - EMERGENTIA Wind Tunnel (wing no sweep)\n10 - EMERGENTIA Wind Tunnel (wing with sweep)\n11 - EMERGENTIA Manufactured\n12 - ALO INTA\n13 - ALO INTA with Fuel Cell\n14 - A320-200\n15 - EMERGENTIA tandem-wing\n16 - TARSIS 75\n17 - TARSIS 120');
dlgtitle = 'Aircraft Design Selection';
dims = [1 35];
definput = {'1'};
opts.Interpreter = 'tex';
answer_CASE_AC = inputdlg(prompt,dlgtitle,dims,definput);
case_AC = str2num(answer_CASE_AC{1});

% Defines path where the MAT files are stored
pathname = fileparts ('ac_models');

switch case_AC
    case 1 % case_AC = 1 - EMERGENTIA 1:1
        OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        OUTPUT_read_XLSX_EMERGENTIA = OUTPUT_read_XLSX;
        save('ac_models/OUTPUT_read_XLSX_EMERGENTIA_100.mat', 'OUTPUT_read_XLSX')
    case 2 % case_AC = 2 - EMERGENTIA Manufactured
        OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        OUTPUT_read_XLSX_EMERGENTIA = OUTPUT_read_XLSX;
%         save OUTPUT_read_XLSX_EMERGENTIA_Mft.mat
        save('ac_models/OUTPUT_read_XLSX_EMERGENTIA_Mft.mat', 'OUTPUT_read_XLSX')
    case 3 % case_AC = 3 - PEPIÑO XXL
        OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        OUTPUT_read_XLSX_PEPINO = OUTPUT_read_XLSX;
%         save OUTPUT_read_XLSX_PEPINO.mat
        save('ac_models/OUTPUT_read_XLSX_PEPINO.mat', 'OUTPUT_read_XLSX')
    case 4 % Commercial example
        OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        OUTPUT_read_XLSX_Comertial = OUTPUT_read_XLSX;
%         save OUTPUT_read_XLSX_Comertial.mat
        save('ac_models/OUTPUT_read_XLSX_A400M.mat', 'OUTPUT_read_XLSX')
    case 5 % WIG Aircraft
        OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        OUTPUT_read_XLSX_WIG = OUTPUT_read_XLSX;
%         save OUTPUT_read_XLSX_WIG.mat
        save('ac_models/OUTPUT_read_XLSX_WIG.mat', 'OUTPUT_read_XLSX')
    case 6 % case_AC = 6 CERVERA
        OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        OUTPUT_read_XLSX_CERVERA = OUTPUT_read_XLSX;
%         save OUTPUT_read_XLSX_CERVERA.mat
        save('ac_models/OUTPUT_read_XLSX_NIMVUS.mat', 'OUTPUT_read_XLSX')
    case 7 % case_AC = 7 - Existing AC - B747
        OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        OUTPUT_read_XLSX_EAC = OUTPUT_read_XLSX;
%         save OUTPUT_read_XLSX_EAC.mat
        save('ac_models/OUTPUT_read_XLSX_EAC.mat', 'OUTPUT_read_XLSX')
    case 8 % case_AC = 8 - E26 Tamiz
        OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        OUTPUT_read_XLSX_Tamiz = OUTPUT_read_XLSX;
%         save OUTPUT_read_XLSX_Tamiz.mat
        save('ac_models/OUTPUT_read_XLSX_Tamiz.mat', 'OUTPUT_read_XLSX')
    case 9 % case_AC = 9 - EMERGENTIA Wind Tunnel Model - (no sweep)
        OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        OUTPUT_read_XLSX_EMERGENTIA_WT2 = OUTPUT_read_XLSX;
%         save OUTPUT_read_XLSX_EMERGENTIA_WT2.mat
        save('ac_models/OUTPUT_read_XLSX_EMERGENTIA_WT2.mat', 'OUTPUT_read_XLSX')
    case 10 % case_AC = 10 - EMERGENTIA Wind Tunnel Model - (sweep)
        OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        OUTPUT_read_XLSX_EMERGENTIA_WT1 = OUTPUT_read_XLSX;
%         save OUTPUT_read_XLSX_EMERGENTIA_WT1.mat
        save('ac_models/OUTPUT_read_XLSX_EMERGENTIA_WT1.mat', 'OUTPUT_read_XLSX')
    case 11 % case_AC = 11 - EMERGENTIA Manufactured
        OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        OUTPUT_read_XLSX_EMERGENTIA_MFT = OUTPUT_read_XLSX;
%         save OUTPUT_read_XLSX_EMERGENTIA_MFT.mat
        save('ac_models/OUTPUT_read_XLSX_EMERGENTIA_MFT.mat', 'OUTPUT_read_XLSX')
    case 12 % case_AC = 12 - ALO
        OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        OUTPUT_read_XLSX_EMERGENTIA_ALO = OUTPUT_read_XLSX;
%         save OUTPUT_read_XLSX_ALO.mat
        save('ac_models/OUTPUT_read_XLSX_ALO.mat', 'OUTPUT_read_XLSX')
    case 13 % case_AC = 13 - ALO Fuel Cell
        OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        OUTPUT_read_XLSX_EMERGENTIA_ALO_FC = OUTPUT_read_XLSX;
%         save OUTPUT_read_XLSX_ALO_FC.mat
        save('ac_models/OUTPUT_read_XLSX_ALO_FC.mat', 'OUTPUT_read_XLSX')
    case 14 % case_AC = 14 - A320-200
        OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        OUTPUT_read_XLSX_EMERGENTIA_A320_FC = OUTPUT_read_XLSX;
%         save OUTPUT_read_XLSX_A320_200.mat
        save('ac_models/OUTPUT_read_XLSX_A320_200.mat', 'OUTPUT_read_XLSX')
    case 15 % case_AC = 15 - EMERGENTIA tandem-wing
        OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        OUTPUT_read_XLSX_EMERGENTIA_TW_FC = OUTPUT_read_XLSX;
%         save OUTPUT_read_XLSX_EMERGENTIA_Tandem_wing.mat
        save('ac_models/OUTPUT_read_XLSX_EMERGENTIA_Tandem_wing.mat', 'OUTPUT_read_XLSX')
    case 16 % case_AC = 15 - Tarsis 75
        OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        OUTPUT_read_XLSX_EMERGENTIA_T95_FC = OUTPUT_read_XLSX;
%         save OUTPUT_read_XLSX_Tarsis75.mat
        save('ac_models/OUTPUT_read_XLSX_Tarsis95.mat', 'OUTPUT_read_XLSX')
    case 17 % case_AC = 15 - Tarsis 120
        OUTPUT_read_XLSX = read_DATA_XLSX_advanced(case_AC);
        OUTPUT_read_XLSX_EMERGENTIA_T120 = OUTPUT_read_XLSX;
%         save OUTPUT_read_XLSX_Tarsis120.mat
        save('ac_models/OUTPUT_read_XLSX_Tarsis120.mat', 'OUTPUT_read_XLSX')
end

if OUTPUT_read_XLSX.STUDY_flags.PERFORMANCE_STUDY == 1

                         % 1st Question
                        prompt4 = sprintf('ENTER MISSION SEGMENTS SELECTION:\n1 - TAXY\n2 - TAKE-OFF\n3 - CLIMB\n4 - VTOL CLIMB\n5 - CRUISE\n6 - LOAD DEPLOYMENT\n7 - TURN\n8 - DESCENT\n9 - VTOL DESCENT\n10 - ALTERNATIVE AIRPORT CLIMB TO 3000FT\n11 - TURN LOITTER 45MIN\n12 - LANDING\n13 - DUMMY\n\nNOTE: Introduce the mission selection as a succession of numbers separated by a blank without any commas or brackets');
                        dlgtitle = 'Mission Segments Selection';
                        dims = [1 70];
                        answer_SEGMENTS = inputdlg(prompt4,dlgtitle,dims); % 3 5 8
                        SEGMENTS = str2num(answer_SEGMENTS{1});
                        SEGMENTS_STYPE = zeros(1,length(SEGMENTS));
                        SEGMENTS_STYPE_InpDat = zeros(length(SEGMENTS),3);
                        flag = 0;
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
                        OUTPUT_read_XLSX.IPP_flags.V_d              = zeros(1,length(SEGMENTS));
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
                        
                        % 2nd Question
                        for i=1:length(SEGMENTS)
                            flag = flag + 1;
                            if SEGMENTS(i) == 1 %TAXY
                                prompt5 = sprintf('ENTER TAXY SEGMENT SUBTYPE SELECTION:\n1 - Giving atmospheric data and taxy configuration');
                                dlgtitle = ['Mission Segment ',num2str(i),' Subtype Selection'];
                                dims = [1 70];
                                answer_SEGMENTS_STYPE = inputdlg(prompt5,dlgtitle,dims);
                                SEGMENTS_STYPE(i) = str2num(answer_SEGMENTS_STYPE{1});
                                % 3rd Question
                                if SEGMENTS_STYPE(i) == 1
                                    prompt6 = sprintf('ENTER TAXY INPUT DATA:\n1 - T [ºC]\n2 - h [m]\n3 - P [Pa]\n4 - delta_T [-]\n5 - V [m/s]\n6 - t [s]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:6) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.temp_local_taxy(flag)  = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_taxy(flag)   = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.P_local_taxy(flag)     = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.delta_T_taxy(flag)     = SEGMENTS_STYPE_InpDat(i,4);
                                    OUTPUT_read_XLSX.IPP_flags.V_taxy(flag)           = SEGMENTS_STYPE_InpDat(i,5);
                                    OUTPUT_read_XLSX.IPP_flags.t_taxy(flag)           = SEGMENTS_STYPE_InpDat(i,6);
                                end
                            end
                            if SEGMENTS(i) == 2 %TAKE-OFF
                                prompt5 = sprintf('ENTER TAKE-OFF SEGMENT SUBTYPE SELECTION:\n1 - Giving atmospheric data and take-off configuration');
                                dlgtitle = ['Mission Segment ',num2str(i),' Subtype Selection'];
                                dims = [1 70];
                                answer_SEGMENTS_STYPE = inputdlg(prompt5,dlgtitle,dims);
                                SEGMENTS_STYPE(i) = str2num(answer_SEGMENTS_STYPE{1});
                                % 3rd Question
                                if SEGMENTS_STYPE(i) == 1
                                    prompt6 = sprintf('ENTER TAKE-OFF INPUT DATA:\n1 - T [ºC]\n2 - h [m]\n3 - P [Pa]\n4 - mu [-]\n5 - h_obstacle [m]\n6 - gamma [º]\n7 - delta_T [-]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:7) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.temp_local_TO(flag)    = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_TO(flag)     = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.P_local_TO(flag)       = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.mu_TO(flag)            = SEGMENTS_STYPE_InpDat(i,4);
                                    OUTPUT_read_XLSX.IPP_flags.h_obstacle_TO(flag)    = SEGMENTS_STYPE_InpDat(i,5);
                                    OUTPUT_read_XLSX.IPP_flags.gamma_climb_TO(flag)   = SEGMENTS_STYPE_InpDat(i,6)*pi/180;
                                    OUTPUT_read_XLSX.IPP_flags.delta_T_TO(flag)       = SEGMENTS_STYPE_InpDat(i,7);
                                end
                            end
                            if SEGMENTS(i) == 3 %CLIMB
                                prompt5 = sprintf('ENTER CLIMB SEGMENT SUBTYPE SELECTION:\n1 - Giving M and gamma\n2 - Giving EAS and gamma\n3 - Giving TAS and gamma\n4 - Giving M and throttle\n5 - Giving EAS and throttle\n6 - Giving TAS and throttle\n7 - Giving Vi, Vf and gamma\n8 - Steppest Climb\n9 - Fastest Climb\n10 - Giving Vi, Vf and throttle');
                                dlgtitle = ['Mission Segment ',num2str(i),' Subtype Selection'];
                                dims = [1 70];
                                answer_SEGMENTS_STYPE = inputdlg(prompt5,dlgtitle,dims);
                                SEGMENTS_STYPE(i) = str2num(answer_SEGMENTS_STYPE{1});
                                % 3rd Question
                                if SEGMENTS_STYPE(i) == 1
                                    prompt6 = sprintf('ENTER CLIMB INPUT DATA:\n1 - Mach Number [-]\n2 - Gamma [deg]\n3 - h_i [m]\n4 - h_f [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:4) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.Mach_cl(flag)      = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.gamma_cl(flag)     = SEGMENTS_STYPE_InpDat(i,2)*pi/180;
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_cl(flag) = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_cl(flag)   = SEGMENTS_STYPE_InpDat(i,4);
                                end
                                if SEGMENTS_STYPE(i) == 2
                                    prompt6 = sprintf('ENTER CLIMB INPUT DATA:\n1 - EAS [m/s]\n2 - Gamma [deg]\n3 - h_i [m]\n4 - h_f [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:4) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.EAS_cl(flag)       = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.gamma_cl(flag)     = SEGMENTS_STYPE_InpDat(i,2)*pi/180;
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_cl(flag) = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_cl(flag)   = SEGMENTS_STYPE_InpDat(i,4);
                                end
                                if SEGMENTS_STYPE(i) == 3
                                    prompt6 = sprintf('ENTER CLIMB INPUT DATA:\n1 - TAS [m/s]\n2 - Gamma [deg]\n3 - h_i [m]\n4 - h_f [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:4) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.TAS_cl(flag)       = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.gamma_cl(flag)     = SEGMENTS_STYPE_InpDat(i,2)*pi/180;
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_cl(flag) = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_cl(flag)   = SEGMENTS_STYPE_InpDat(i,4);
                                end
                                if SEGMENTS_STYPE(i) == 4
                                    prompt6 = sprintf('ENTER CLIMB INPUT DATA:\n1 - Mach Number [-]\n2 - Throttle [-]\n3 - h_i [m]\n4 - h_f [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:4) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.Mach_cl(flag)      = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.delta_T_cl(flag)   = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_cl(flag) = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_cl(flag)   = SEGMENTS_STYPE_InpDat(i,4);
                                end
                                if SEGMENTS_STYPE(i) == 5
                                    prompt6 = sprintf('ENTER CLIMB INPUT DATA:\n1 - EAS [m/s]\n2 - Throttle [-]\n3 - h_i [m]\n4 - h_f [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:4) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.EAS_cl(flag)       = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.delta_T_cl(flag)   = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_cl(flag) = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_cl(flag)   = SEGMENTS_STYPE_InpDat(i,4);
                                end
                                if SEGMENTS_STYPE(i) == 6
                                    prompt6 = sprintf('ENTER CLIMB INPUT DATA:\n1 - TAS [m/s]\n2 - Throttle [-]\n3 - h_i [m]\n4 - h_f [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:4) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.TAS_cl(flag)       = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.delta_T_cl(flag)   = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_cl(flag) = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_cl(flag)   = SEGMENTS_STYPE_InpDat(i,4);
                                end
                                if SEGMENTS_STYPE(i) == 7
                                    prompt6 = sprintf('ENTER CLIMB INPUT DATA:\n1 - Vi [m/s]\n2 - Vf [m/s]\n3 - Gamma [deg]\n4 - h_i [m]\n5 - h_f [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:5) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.V_ini_cl(flag)     = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.V_fin_cl(flag)     = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.gamma_cl(flag)     = SEGMENTS_STYPE_InpDat(i,3)*pi/180;
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_cl(flag) = SEGMENTS_STYPE_InpDat(i,4);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_cl(flag)   = SEGMENTS_STYPE_InpDat(i,5);
                                end
                                if SEGMENTS_STYPE(i) == 8
                                    prompt6 = sprintf('ENTER CLIMB INPUT DATA:\n1 - Throttle [-]\n2 - h_i [m]\n3 - h_f [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:3) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.delta_T_cl(flag)     = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_cl(flag)   = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_cl(flag)     = SEGMENTS_STYPE_InpDat(i,3);
                                end
                                if SEGMENTS_STYPE(i) == 9
                                    prompt6 = sprintf('ENTER CLIMB INPUT DATA:\n1 - Throttle [-]\n2 - h_i [m]\n3 - h_f [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:3) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.delta_T_cl(flag)     = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_cl(flag)   = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_cl(flag)     = SEGMENTS_STYPE_InpDat(i,3);
                                end
                                if SEGMENTS_STYPE(i) == 10
                                    prompt6 = sprintf('ENTER CLIMB INPUT DATA:\n1 - Vi [m/s]\n2 - Vf [m/s]\n3 - Throttle [-]\n4 - h_i [m]\n5 - h_f [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:5) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.V_ini_cl(flag)     = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.V_fin_cl(flag)     = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.delta_T_cl(flag)   = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_cl(flag) = SEGMENTS_STYPE_InpDat(i,4);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_cl(flag)   = SEGMENTS_STYPE_InpDat(i,5);
                                end
                            end
                            if SEGMENTS(i) == 4 %VTOL CLIMB
                                prompt5 = sprintf('ENTER VTOL CLIMB SEGMENT SUBTYPE SELECTION:\n1 - Giving hf and throttle\n2 - Hovering giving t\n3 - Giving hf and vclimb\n4 - Hovering giving mbat');
                                dlgtitle = ['Mission Segment ',num2str(i),' Subtype Selection'];
                                dims = [1 70];
                                answer_SEGMENTS_STYPE = inputdlg(prompt5,dlgtitle,dims);
                                SEGMENTS_STYPE(i) = str2num(answer_SEGMENTS_STYPE{1});
                                % 3rd Question
                                if SEGMENTS_STYPE(i) == 1
                                    prompt6 = sprintf('ENTER VTOL CLIMB INPUT DATA:\n1 - h_i [m]\n2 - h_f [m]\n3 - Throttle [-]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:3) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_vtcl(flag) = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_vtcl(flag)   = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.delta_T_vtcl(flag)   = SEGMENTS_STYPE_InpDat(i,3);
                                end
                                if SEGMENTS_STYPE(i) == 2
                                    prompt6 = sprintf('ENTER VTOL CLIMB INPUT DATA:\n1 - h_i [m]\n2 - t_hover [s]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:2) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_vtcl(flag) = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.t_hover(flag)        = SEGMENTS_STYPE_InpDat(i,2);
                                end
                                if SEGMENTS_STYPE(i) == 3
                                    prompt6 = sprintf('ENTER VTOL CLIMB INPUT DATA:\n1 - h_i [m]\n2 - h_f [m]\n3 - vclimb [m/s]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:3) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_vtcl(flag) = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_vtcl(flag)   = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.vclimb_vtcl(flag)    = SEGMENTS_STYPE_InpDat(i,3);
                                end
                                if SEGMENTS_STYPE(i) == 4
                                    prompt6 = sprintf('ENTER VTOL CLIMB INPUT DATA:\n1 - h_i [m]\n2 - mbat [kg]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:2) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_vtcl(flag) = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.mbat_vtcl(flag)      = SEGMENTS_STYPE_InpDat(i,2);
                                end
                            end
                            if SEGMENTS(i) == 5 %CRUISE
                                prompt5 = sprintf('ENTER CRUISE SEGMENT SUBTYPE SELECTION:\n1 - Giving M and distance\n2 - Giving CL and distance\n3 - Giving Vi, Vf and throttle\n4 - Giving M and CD(M)\n5 - Max range giving Wf\n6 - Max endurance giving Wf\n7 - Giving M and m_bat');
                                dlgtitle = ['Mission Segment ',num2str(i),' Subtype Selection'];
                                dims = [1 70];
                                answer_SEGMENTS_STYPE = inputdlg(prompt5,dlgtitle,dims);
                                SEGMENTS_STYPE(i) = str2num(answer_SEGMENTS_STYPE{1});
                                % 3rd Question
                                if SEGMENTS_STYPE(i) == 1
                                    prompt6 = sprintf('ENTER CRUISE INPUT DATA:\n1 - V_cr [m/s]\n2 - Distance [m]\n3 - h_i [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:3) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.V_cr(flag)           = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.dist_final_cr(flag)  = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_cr(flag)   = SEGMENTS_STYPE_InpDat(i,3);
                                end
                                if SEGMENTS_STYPE(i) == 2
                                    prompt6 = sprintf('ENTER CRUISE INPUT DATA:\n1 - Distance [m]\n2 - h_i [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:2) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.dist_final_cr(flag)  = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_cr(flag)   = SEGMENTS_STYPE_InpDat(i,2);
                                end
                                if SEGMENTS_STYPE(i) == 3
                                    prompt6 = sprintf('ENTER CRUISE INPUT DATA:\n1 - Vi [m/s]\n2 - Vf [m/s]\n3 - Throttle [-]\n4 - Distance [m]\n5 - h_i [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:5) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.V_ini_cr(flag)           = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.V_fin_cr(flag)           = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.delta_T_cr(flag)         = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.dist_final_cr(flag)      = SEGMENTS_STYPE_InpDat(i,4);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_cr(flag)       = SEGMENTS_STYPE_InpDat(i,5);
                                end
                                if SEGMENTS_STYPE(i) == 4
                                    prompt6 = sprintf('ENTER CRUISE INPUT DATA:\n1 - V_cr [m/s]\n2 - CD_0 [-]\n3 - K1 [-]\n4 - K2 [-]\n5 - Distance [m]\n6 - h_i [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:6) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.V_cr(flag)               = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.Cd0_cr(flag)             = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.k1_cr(flag)              = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.k2_cr(flag)              = SEGMENTS_STYPE_InpDat(i,4);
                                    OUTPUT_read_XLSX.IPP_flags.dist_final_cr(flag)      = SEGMENTS_STYPE_InpDat(i,5);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_cr(flag)       = SEGMENTS_STYPE_InpDat(i,6);
                                end
                                if SEGMENTS_STYPE(i) == 5
                                    prompt6 = sprintf('ENTER CRUISE INPUT DATA:\n1 - Wf [kg]\n2 - h_i [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:2) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.fuel_cr(flag)        = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_cr(flag)   = SEGMENTS_STYPE_InpDat(i,2);
                                end
                                if SEGMENTS_STYPE(i) == 6
                                    prompt6 = sprintf('ENTER CRUISE INPUT DATA:\n1 - Wf [kg]\n2 - h_i [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:2) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.fuel_cr(flag)        = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_cr(flag)   = SEGMENTS_STYPE_InpDat(i,2);
                                end
                                if SEGMENTS_STYPE(i) == 7
                                    prompt6 = sprintf('ENTER CRUISE INPUT DATA:\n1 - V_cr [m/s]\n2 - m_bat [kg]\n3 - h_i [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:3) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.V_cr(flag)           = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.mbat(flag)           = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_cr(flag)   = SEGMENTS_STYPE_InpDat(i,3);
                                end
                            end
                            if SEGMENTS(i) == 6 %LOAD DEPLOYMENT
                                prompt5 = sprintf('ENTER LOAD DEPLOYMENT SEGMENT SUBTYPE SELECTION:\n1 - Giving the value of the load');
                                dlgtitle = ['Mission Segment ',num2str(i),' Subtype Selection'];
                                dims = [1 70];
                                answer_SEGMENTS_STYPE = inputdlg(prompt5,dlgtitle,dims);
                                SEGMENTS_STYPE(i) = str2num(answer_SEGMENTS_STYPE{1});
                                % 3rd Question
                                if SEGMENTS_STYPE(i) == 1
                                    prompt6 = sprintf('ENTER LOAD DEPLOYMENT INPUT DATA:\n1 - L [kg]');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.carga_loadep(flag) = SEGMENTS_STYPE_InpDat(i);
                                end
                            end
                            if SEGMENTS(i) == 7 %TURN
                                prompt5 = sprintf('ENTER TURN SEGMENT SUBTYPE SELECTION:\n1 - Giving V and throttle\n2 - Giving V and CL\n3 - Giving V and roll angle\n4 - Giving V and load factor\n5 - Giving V and turn radius\n6 - Giving V and yaw rate\n7 - Giving throttle with maximum load factor\n8 - Giving throttle with maximum yaw rate\n9 - Giving throttle with minimum turn radius');
                                dlgtitle = ['Mission Segment ',num2str(i),' Subtype Selection'];
                                dims = [1 70];
                                answer_SEGMENTS_STYPE = inputdlg(prompt5,dlgtitle,dims);
                                SEGMENTS_STYPE(i) = str2num(answer_SEGMENTS_STYPE{1});
                                % 3rd Question
                                if SEGMENTS_STYPE(i) == 1
                                    prompt6 = sprintf('ENTER TURN INPUT DATA:\n1 - V_tr [m/s]\n2 - Throttle [-]\n3 - h_i [m]\n4 - t [s]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:4) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.V_turn(flag)         = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.delta_T_tr(flag)     = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_tr(flag)   = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.t_final_tr(flag)     = SEGMENTS_STYPE_InpDat(i,4);
                                end
                                if SEGMENTS_STYPE(i) == 2
                                    prompt6 = sprintf('ENTER TURN INPUT DATA:\n1 - V_tr [m/s]\n2 - h_i [m]\n3 - t [s]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,3) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.V_turn(flag)         = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_tr(flag)   = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.t_final_tr(flag)     = SEGMENTS_STYPE_InpDat(i,3);
                                end
                                if SEGMENTS_STYPE(i) == 3
                                    prompt6 = sprintf('ENTER TURN INPUT DATA:\n1 - V_tr [m/s]\n2 - Phi [deg]\n3 - h_i [m]\n4 - t [s]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:4) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.V_turn(flag)         = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.phi_tr(flag)         = SEGMENTS_STYPE_InpDat(i,2)*pi/180;
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_tr(flag)   = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.t_final_tr(flag)     = SEGMENTS_STYPE_InpDat(i,4);
                                end
                                if SEGMENTS_STYPE(i) == 4
                                    prompt6 = sprintf('ENTER TURN INPUT DATA:\n1 - V_tr [m/s]\n2 - n [-]\n3 - h_i [m]\n4 - t [s]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:4) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.V_turn(flag)         = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.n_tr(flag)           = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_tr(flag)   = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.t_final_tr(flag)     = SEGMENTS_STYPE_InpDat(i,4);
                                end
                                if SEGMENTS_STYPE(i) == 5
                                    prompt6 = sprintf('ENTER TURN INPUT DATA:\n1 - V_tr [m/s]\n2 - R_tr [m]\n3 - h_i [m]\n4 - t [s]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:4) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.V_turn(flag)         = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.R_tr(flag)           = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_tr(flag)   = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.t_final_tr(flag)     = SEGMENTS_STYPE_InpDat(i,4);
                                end
                                if SEGMENTS_STYPE(i) == 6
                                    prompt6 = sprintf('ENTER TURN INPUT DATA:\n1 - V_tr [m/s]\n2 - V_psi [deg/s]\n3 - h_i [m]\n4 - t [s]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:4) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.V_turn(flag)         = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.V_psi(flag)          = SEGMENTS_STYPE_InpDat(i,2)*pi/180;
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_tr(flag)   = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.t_final_tr(flag)     = SEGMENTS_STYPE_InpDat(i,4);
                                end
                                if SEGMENTS_STYPE(i) == 7
                                    prompt6 = sprintf('ENTER TURN INPUT DATA:\n1 - Throttle [-]\n2 - h_i [m]\n3 - t [s]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:3) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.delta_T_tr(flag)     = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_tr(flag)   = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.t_final_tr(flag)     = SEGMENTS_STYPE_InpDat(i,3);
                                end
                                if SEGMENTS_STYPE(i) == 8
                                    prompt6 = sprintf('ENTER TURN INPUT DATA:\n1 - Throttle [-]\n2 - h_i [m]\n3 - t [s]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:3) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.delta_T_tr(flag)     = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_tr(flag)   = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.t_final_tr(flag)     = SEGMENTS_STYPE_InpDat(i,3);
                                end
                                if SEGMENTS_STYPE(i) == 9
                                    prompt6 = sprintf('ENTER TURN INPUT DATA:\n1 - Throttle [-]\n2 - h_i [m]\n3 - t [s]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:3) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.delta_T_tr(flag)     = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_tr(flag)   = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.t_final_tr(flag)     = SEGMENTS_STYPE_InpDat(i,3);
                                end
                            end
                            if SEGMENTS(i) == 8 %DESCENT
                                prompt5 = sprintf('ENTER DESCENT SEGMENT SUBTYPE SELECTION:\n1 - Giving M and gamma\n2 - Giving EAS and gamma\n3 - Giving TAS and gamma\n4 - Giving M and throttle\n5 - Giving EAS and throttle\n6 - Giving TAS and throttle\n7 - Giving Vi, Vf and gamma\n8 - Minimum gamma\n9 - Slowest sink\n10 - Giving Vi, Vf and throttle');
                                dlgtitle = ['Mission Segment ',num2str(i),' Subtype Selection'];
                                dims = [1 70];
                                answer_SEGMENTS_STYPE = inputdlg(prompt5,dlgtitle,dims);
                                SEGMENTS_STYPE(i) = str2num(answer_SEGMENTS_STYPE{1});
                                % 3rd Question
                                if SEGMENTS_STYPE(i) == 1
                                    prompt6 = sprintf('ENTER DESCENT INPUT DATA:\n1 - V_d [m/s]\n2 - Gamma [deg]\n3 - h_i [m]\n4 - h_f [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:4) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.V_d(flag)         = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.gamma_d(flag)     = SEGMENTS_STYPE_InpDat(i,2)*pi/180;
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_d(flag) = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_d(flag)   = SEGMENTS_STYPE_InpDat(i,4);
                                end
                                if SEGMENTS_STYPE(i) == 2
                                    prompt6 = sprintf('ENTER DESCENT INPUT DATA:\n1 - EAS [m/s]\n2 - Gamma [-]\n3 - h_i [m]\n4 - h_f [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:4) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.EAS_d(flag)       = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.gamma_d(flag)     = SEGMENTS_STYPE_InpDat(i,2)*pi/180;
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_d(flag) = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_d(flag)   = SEGMENTS_STYPE_InpDat(i,4);
                                end
                                if SEGMENTS_STYPE(i) == 3
                                    prompt6 = sprintf('ENTER DESCENT INPUT DATA:\n1 - TAS [m/s]\n2 - Gamma [deg]\n3 - h_i [m]\n4 - h_f [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:4) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.TAS_d(flag)       = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.gamma_d(flag)     = SEGMENTS_STYPE_InpDat(i,2)*pi/180;
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_d(flag) = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_d(flag)   = SEGMENTS_STYPE_InpDat(i,4);
                                end
                                if SEGMENTS_STYPE(i) == 4
                                    prompt6 = sprintf('ENTER DESCENT INPUT DATA:\n1 - V_d [m/s]\n2 - Throttle [-]\n3 - h_i [m]\n4 - h_f [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:4) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.V_d(flag)         = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.delta_T_d(flag)   = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_d(flag) = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_d(flag)   = SEGMENTS_STYPE_InpDat(i,4);
                                end
                                if SEGMENTS_STYPE(i) == 5
                                    prompt6 = sprintf('ENTER DESCENT INPUT DATA:\n1 - EAS [m/s]\n2 - Throttle [-]\n3 - h_i [m]\n4 - h_f [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:4) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.EAS_d(flag)       = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.delta_T_d(flag)   = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_d(flag) = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_d(flag)   = SEGMENTS_STYPE_InpDat(i,4);
                                end
                                if SEGMENTS_STYPE(i) == 6
                                    prompt6 = sprintf('ENTER DESCENT INPUT DATA:\n1 - TAS [m/s]\n2 - Throttle [-]\n3 - h_i [m]\n4 - h_f [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:4) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.TAS_d(flag)       = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.delta_T_d(flag)   = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_d(flag) = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_d(flag)   = SEGMENTS_STYPE_InpDat(i,4);
                                end
                                if SEGMENTS_STYPE(i) == 7
                                    prompt6 = sprintf('ENTER DESCENT INPUT DATA:\n1 - Vi [m/s]\n2 - Vf [m/s]\n3 - Gamma [deg]\n4 - h_i [m]\n5 - h_f [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:5) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.V_ini_d(flag)     = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.V_fin_d(flag)     = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.gamma_d(flag)     = SEGMENTS_STYPE_InpDat(i,3)*pi/180;
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_d(flag) = SEGMENTS_STYPE_InpDat(i,4);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_d(flag)   = SEGMENTS_STYPE_InpDat(i,5);
                                end
                                if SEGMENTS_STYPE(i) == 8
                                    prompt6 = sprintf('ENTER DESCENT INPUT DATA:\n1 - Throttle [-]\n2 - h_i [m]\n3 - h_f [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:3) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.delta_T_d(flag)     = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_d(flag)   = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_d(flag)     = SEGMENTS_STYPE_InpDat(i,3);
                                end
                                if SEGMENTS_STYPE(i) == 9
                                    prompt6 = sprintf('ENTER DESCENT INPUT DATA:\n1 - Throttle [-]\n2 - h_i [m]\n3 - h_f [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:3) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.delta_T_d(flag)     = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_d(flag)   = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_d(flag)     = SEGMENTS_STYPE_InpDat(i,3);
                                end                               
                                if SEGMENTS_STYPE(i) == 10
                                    prompt6 = sprintf('ENTER DESCENT INPUT DATA:\n1 - Vi [m/s]\n2 - Vf [m/s]\n3 - Throttle [-]\n4 - h_i [m]\n5 - h_f [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:5) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.V_ini_d(flag)     = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.V_fin_d(flag)     = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.delta_T_d(flag)   = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_d(flag) = SEGMENTS_STYPE_InpDat(i,4);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_d(flag)   = SEGMENTS_STYPE_InpDat(i,5);
                                end
                            end
                            if SEGMENTS(i) == 9 %VTOL DESCENT
                                prompt5 = sprintf('ENTER VTOL DESCENT SEGMENT SUBTYPE SELECTION:\n1 - Giving hi, hf and vdes');
                                dlgtitle = ['Mission Segment ',num2str(i),' Subtype Selection'];
                                dims = [1 70];
                                answer_SEGMENTS_STYPE = inputdlg(prompt5,dlgtitle,dims);
                                SEGMENTS_STYPE(i) = str2num(answer_SEGMENTS_STYPE{1});
                                % 3rd Question
                                if SEGMENTS_STYPE(i) == 1
                                    prompt6 = sprintf('ENTER VTOL DESCENT INPUT DATA:\n1 - h_i [m]\n2 - h_f [m]\n3 - v_des [m/s]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:3) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_vtd(flag) = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_vtd(flag)   = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.vdes_vtd(flag)     = SEGMENTS_STYPE_InpDat(i,3);
                                end
                            end
                            if SEGMENTS(i) == 10 %ALTERNATE AIRPORT CLIMB (3000 ft)
                                prompt5 = sprintf('ENTER CLIMB SEGMENT SUBTYPE SELECTION:\n1 - Giving M and gamma\n2 - Giving EAS and gamma\n3 - Giving TAS and gamma\n4 - Giving M and throttle\n5 - Giving EAS and throttle\n6 - Giving TAS and throttle\n7 - Giving Vi, Vf and gamma\n8 - Steppest Climb\n9 - Fastest Climb\n10 - Giving Vi, Vf and throttle');
                                dlgtitle = ['Mission Segment ',num2str(i),' Subtype Selection'];
                                dims = [1 70];
                                answer_SEGMENTS_STYPE = inputdlg(prompt5,dlgtitle,dims);
                                SEGMENTS_STYPE(i) = str2num(answer_SEGMENTS_STYPE{1});
                                % 3rd Question
                                if SEGMENTS_STYPE(i) == 1
                                    prompt6 = sprintf('ENTER CLIMB INPUT DATA:\n1 - Mach Number [-]\n2 - Gamma [deg]\n3 - h_i [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:3) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.Mach_cl(flag)      = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.gamma_cl(flag)     = SEGMENTS_STYPE_InpDat(i,2)*pi/180;
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_cl(flag) = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_cl(flag)   = 914.4; %3000 ft
                                end
                                if SEGMENTS_STYPE(i) == 2
                                    prompt6 = sprintf('ENTER CLIMB INPUT DATA:\n1 - EAS [m/s]\n2 - Gamma [deg]\n3 - h_i [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:3) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.EAS_cl(flag)       = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.gamma_cl(flag)     = SEGMENTS_STYPE_InpDat(i,2)*pi/180;
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_cl(flag) = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_cl(flag)   = 914.4; %3000 ft
                                end
                                if SEGMENTS_STYPE(i) == 3
                                    prompt6 = sprintf('ENTER CLIMB INPUT DATA:\n1 - TAS [m/s]\n2 - Gamma [deg]\n3 - h_i [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:3) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.TAS_cl(flag)       = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.gamma_cl(flag)     = SEGMENTS_STYPE_InpDat(i,2)*pi/180;
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_cl(flag) = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_cl(flag)   = 914.4; %3000 ft
                                end
                                if SEGMENTS_STYPE(i) == 4
                                    prompt6 = sprintf('ENTER CLIMB INPUT DATA:\n1 - Mach Number [-]\n2 - Throttle [-]\n3 - h_i [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:3) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.Mach_cl(flag)      = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.delta_T_cl(flag)   = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_cl(flag) = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_cl(flag)   = 914.4; %3000 ft
                                end
                                if SEGMENTS_STYPE(i) == 5
                                    prompt6 = sprintf('ENTER CLIMB INPUT DATA:\n1 - EAS [m/s]\n2 - Throttle [-]\n3 - h_i [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:4) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.EAS_cl(flag)       = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.delta_T_cl(flag)   = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_cl(flag) = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_cl(flag)   = 914.4; %3000 ft
                                end
                                if SEGMENTS_STYPE(i) == 6
                                    prompt6 = sprintf('ENTER CLIMB INPUT DATA:\n1 - TAS [m/s]\n2 - Throttle [-]\n3 - h_i [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:4) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.TAS_cl(flag)       = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.delta_T_cl(flag)   = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_cl(flag) = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_cl(flag)   = 914.4; %3000 ft
                                end
                                if SEGMENTS_STYPE(i) == 7
                                    prompt6 = sprintf('ENTER CLIMB INPUT DATA:\n1 - Vi [m/s]\n2 - Vf [m/s]\n3 - Gamma [deg]\n4 - h_i [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:4) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.V_ini_cl(flag)     = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.V_fin_cl(flag)     = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.gamma_cl(flag)     = SEGMENTS_STYPE_InpDat(i,3)*pi/180;
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_cl(flag) = SEGMENTS_STYPE_InpDat(i,4);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_cl(flag)   = 914.4; %3000 ft
                                end
                                if SEGMENTS_STYPE(i) == 8
                                    prompt6 = sprintf('ENTER CLIMB INPUT DATA:\n1 - Throttle [-]\n2 - h_i [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:2) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.delta_T_cl(flag)     = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_cl(flag)   = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_cl(flag)     = 914.4; %3000 ft
                                end
                                if SEGMENTS_STYPE(i) == 9
                                    prompt6 = sprintf('ENTER CLIMB INPUT DATA:\n1 - Throttle [-]\n2 - h_i [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:2) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.delta_T_cl(flag)     = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_cl(flag)   = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_cl(flag)     = 914.4; %3000 ft
                                end
                                if SEGMENTS_STYPE(i) == 10
                                    prompt6 = sprintf('ENTER CLIMB INPUT DATA:\n1 - Vi [m/s]\n2 - Vf [m/s]\n3 - Throttle [-]\n4 - h_i [m]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:4) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.V_ini_cl(flag)     = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.V_fin_cl(flag)     = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.delta_T_cl(flag)   = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_cl(flag) = SEGMENTS_STYPE_InpDat(i,4);
                                    OUTPUT_read_XLSX.IPP_flags.h_final_cl(flag)   = 914.4; %3000 ft
                                end
                            end
                            if SEGMENTS(i) == 12 %LANDING
                                prompt5 = sprintf('ENTER LANDING SEGMENT SUBTYPE SELECTION:\n1 - Giving atmospheric data and landing configuration');
                                dlgtitle = ['Mission Segment ',num2str(i),' Subtype Selection'];
                                dims = [1 70];
                                answer_SEGMENTS_STYPE = inputdlg(prompt5,dlgtitle,dims);
                                SEGMENTS_STYPE(i) = str2num(answer_SEGMENTS_STYPE{1});
                                % 3rd Question
                                if SEGMENTS_STYPE(i) == 1
                                    prompt6 = sprintf('ENTER LANDING INPUT DATA:\n1 - T [ºC]\n2 - h [m]\n3 - P [Pa]\n4 - mu [-]\n5 - delta_T [-]\n6 - t_brake [s]\n\nNOTE: Introduce the input data as a succession of numbers separated by a blank without any commas or brackets');
                                    dlgtitle = ['Mission Segment ',num2str(i),' Input Data'];
                                    dims = [1 70];
                                    answer_SEGMENT_InputData = inputdlg(prompt6,dlgtitle,dims);
                                    SEGMENTS_STYPE_InpDat(i,1:6) = str2num(answer_SEGMENT_InputData{1});
                                    OUTPUT_read_XLSX.IPP_flags.temp_local_LND(flag)   = SEGMENTS_STYPE_InpDat(i,1);
                                    OUTPUT_read_XLSX.IPP_flags.h_inicial_LND(flag)    = SEGMENTS_STYPE_InpDat(i,2);
                                    OUTPUT_read_XLSX.IPP_flags.P_local_LND(flag)      = SEGMENTS_STYPE_InpDat(i,3);
                                    OUTPUT_read_XLSX.IPP_flags.mu_LND(flag)           = SEGMENTS_STYPE_InpDat(i,4);
                                    OUTPUT_read_XLSX.IPP_flags.delta_T_LND(flag)      = SEGMENTS_STYPE_InpDat(i,5);
                                    OUTPUT_read_XLSX.IPP_flags.t_brake(flag)          = SEGMENTS_STYPE_InpDat(i,6);
                                end
                            end
                        end
                        % Saving Mission Configuration
                        OUTPUT_read_XLSX.PerforMisionSelection_flags.type_missions_WF       = SEGMENTS;
                        OUTPUT_read_XLSX.PerforMisionSelection_flags.num_missions_WF        = length(SEGMENTS);
                        OUTPUT_read_XLSX.PerforMisionSelection_flags.sub_type_missions_WF   = SEGMENTS_STYPE;
%                         OUTPUT_read_XLSX.PerforMisionSelection_flags.climb_mode         = SEGMENTS_STYPE(1);
%                         OUTPUT_read_XLSX.PerforMisionSelection_flags.cruise_mode        = SEGMENTS_STYPE(2);
%                         OUTPUT_read_XLSX.PerforMisionSelection_flags.descent_mode       = SEGMENTS_STYPE(3);
                        % Saving .mat
end
