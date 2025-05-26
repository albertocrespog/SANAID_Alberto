function DeltaDrag = get_DeltaDrag_Asymmetry(case_AC,conditions,OUTPUT_read_XLSX,Storing_AERO_DATA,Geo_tier,rho,V)
 
% Geometric data
S_w1 = Geo_tier.S_w1;
b_w1 = Geo_tier.b_w1;
q_inf = 0.5*rho*V^2;

switch case_AC
    case 1 % case_AC = 1 - EMERGENTIA 1:1
        D_wing1 = 0; % Left Wing
        D_wing2 = 0; % Right Wing
    case 2 % case_AC = 2 - EMERGENTIA 1:2
        D_wing1 = 0; % Left Wing
        D_wing2 = 0; % Right Wing
    case 3 % case_AC = 3 - PEPIÑO XXL
        D_wing1 = 0; % Left Wing
        D_wing2 = 0; % Right Wing
    case 4 % case_AC = 4 - COMERCIAL
        D_wing1 = 0; % Left Wing
        D_wing2 = 0; % Right Wing
    case 5 % case_AC = 5 - WIGL
        D_wing1 = 0; % Left Wing
        D_wing2 = 0; % Right Wing
    case 6 % case_AC = 6 - CERVERA
        D_wing1 = 0; % Left Wing
        D_wing2 = 0; % Right Wing
    case 7 % Existing Aircraft = 7
        D_wing1 = 0; % Left Wing
        D_wing2 = 0; % Right Wing
    case 8 % case_AC = 8 - TAMIZ
        D_wing1 = 0; % Left Wing
        D_wing2 = 0; % Right Wing
    case 9 % case_AC = 9 - EMERGENTIA Wind Tunnel Model - w2
        D_wing1 = 0; % Left Wing
        D_wing2 = 0; % Right Wing
    case 10 % case_AC = 10 - EMERGENTIA Wind Tunnel Model - w1 (sweep)
        D_wing1 = 0; % Left Wing
        D_wing2 = 0; % Right Wing
    case 11 % case_AC = 11 - EMERGENTIA Manufacturing
        D_wing1 = 0; % Left Wing
        D_wing2 = 0; % Right Wing
    case 12 % case_AC = 12 - ALO
        D_wing1 = 0; % Left Wing
        D_wing2 = 0; % Right Wing
    case 13 % case_AC = 13 - ALO Fuel Cell
        D_wing1 = 0; % Left Wing
        D_wing2 = 0; % Right Wing
    case 14 % case_AC = 14 - A320
        D_wing1 = 0; % Left Wing
        D_wing2 = 0; % Right Wing
    case 15 % case_AC = 15 - EMERGENTIA Tandem Wing
        D_wing1 = 0; % Left Wing
        D_wing2 = 0; % Right Wing
    case 16 % case_AC = 16 - Tarsis T75
        D_wing1 = 0; % Left Wing
        D_wing2 = 0; % Right Wing
    case 17 % case_AC = 17 - Tarsis T120
        % Tarsis 120
        if conditions.T120 == 1
%             N_moment = conditions.conditions_TRIM_lat.N_moment;
            % Drag asymmetry
            CD0_missile = Storing_AERO_DATA.Aero_TH.CD0_missile;
            CD0_pod = Storing_AERO_DATA.Aero_TH.CD0_pod;
            D_missile = q_inf*S_w1*CD0_missile;
            D_pod = q_inf*S_w1*CD0_pod;           
            % Determines the drag that generates Fox configuration
            switch conditions.n_MSL
                % The missiles are selected such that the asymmetry in
                % Drag produces a positive yawing moment
                case 0
                    DeltaDrag_REF = 0*D_missile;
                    D_wing1 = DeltaDrag_REF; % Left Wing
                    D_wing2 = DeltaDrag_REF; % Right Wing
                case 1
%                     DeltaDrag = 1*D_missile + D_pod; NO POD?????
                    DeltaDrag_REF = 1*D_missile;
                    D_wing1 = 0; % Left Wing
                    D_wing2 = DeltaDrag_REF; % Right Wing
                case 2
%                     DeltaDrag = 2*D_missile + D_pod; NO POD?????
                    DeltaDrag_REF = 2*D_missile;
                    D_wing1 = 0; % Left Wing
                    D_wing2 = DeltaDrag_REF; % Right Wing

                case 3
%                     DeltaDrag = 1*D_missile + D_pod; NO POD?????
                    DeltaDrag_REF = 1*D_missile;
                    D_wing1 = 0; % Left Wing
                    D_wing2 = DeltaDrag_REF; % Right Wing
                case 4
                    DeltaDrag_REF = 0*D_missile;
                    D_wing1 = DeltaDrag_REF; % Left Wing
                    D_wing2 = DeltaDrag_REF; % Right Wing
            end
        
        elseif conditions.T120 == 2 || conditions.T120 == 3 || conditions.T120 == 4%% Tarsis 120 + KSA
%             N_moment = conditions.conditions_TRIM_lat.N_moment;
            % Drag asymmetry
            CD0_missile = Storing_AERO_DATA.Aero_TH.CD0_missile;
            CD0_pod = Storing_AERO_DATA.Aero_TH.CD0_pod;
            D_missile = q_inf*S_w1*CD0_missile;
            D_pod = q_inf*S_w1*CD0_pod;
            switch conditions.n_MSL
                % The missiles are selected such that the asymmetry in
                % Drag produces a positive yawing moment. assume that KSA
                % is placed in right wing
               case 0
                    DeltaDrag_REF = 0*D_missile;
                    D_wing1 = DeltaDrag_REF; % Left Wing
                    D_wing2 = DeltaDrag_REF; % Right Wing
                case 3 % Solo Bola
%                     DeltaDrag = 1*D_missile + D_pod; NO POD?????
                    DeltaDrag_REF = 1*D_missile;
                    D_wing1 = 0; % Left Wing
                    D_wing2 = DeltaDrag_REF; % Right Wing
                case 4 % KSA Completo BAT + Bola
                    DeltaDrag_REF = 1*D_missile;
                    D_wing1 = DeltaDrag_REF; % Left Wing
                    D_wing2 = DeltaDrag_REF; % Right Wing
            end
        end
%         DeltaDrag = 2*D_missile + D_pod;
%         DeltaNDrag = DeltaDrag*N_moment; % Yawing Moment due to drag asymetry
%         DeltaNDrag = -D_wing1*N_moment + D_wing2*N_moment; % Yawing Moment due to drag asymetry
    case 18 % case_AC = 18 - BAT
        D_wing1 = 0; % Left Wing
        D_wing2 = 0; % Right Wing
    case 19 % case_AC = 19 - EVOL1
        D_wing1 = 0; % Left Wing
        D_wing2 = 0; % Right Wing
    case 20 % case_AC = 20 - EVOL2
        D_wing1 = 0; % Left Wing
        D_wing2 = 0; % Right Wing
    case 21 % case_AC = 21 - EVOL3
        D_wing1 = 0; % Left Wing
        D_wing2 = 0; % Right Wing
end
DeltaDrag.D_wing1 = D_wing1;
DeltaDrag.D_wing2 = D_wing2;
% DeltaDrag.DeltaNDrag = DeltaNDrag;