function [X,Y,Z_study_h,Fig] = Sensitivity_Study_Horizontal_Flight(Aero_TH,SF_prop,RPMMAX_APC,Prop_data,Weight_tier,conv_UNITS,...
    Geo_tier,AC_CONFIGURATION,Fig,PLOTS_SENSITIVITY_HORIZONTAL,Performance,solve_w_fero,...
    N_prop,pp_D_prop_min,pp_D_prop_max,N_contour_lines,PLOTS_Prop_Performance,Prop_selection,prefix,Plot_Options)

i = 0;
D_prop = Prop_data.D_prop;
S_ref = Geo_tier.S_ref;
rho = Performance.rho;
g = conv_UNITS.g;
m_TOW = Weight_tier.m_TOW;
W_MTOW = m_TOW*g;

SAVE_FIGS = Plot_Options.SAVE_FIGS;
MATLAB_in = Plot_Options.MATLAB_in;

n_eng = AC_CONFIGURATION.n_eng;  % Número de motores
D_prop_vect=linspace(pp_D_prop_min*D_prop,pp_D_prop_max*D_prop,N_prop);
dist = Performance.Range;    % Distancia de la misión
time = Performance.Endurance*60;

type_battery = AC_CONFIGURATION.type_battery;
switch type_battery
    case 1
        SE_LiFePO4 = Weight_tier.M_ENERGY.SE_LiFePO4; % Wh/kg Specific Energy  90–160 Wh/kg (320–580 J/g or kJ/kg)
        SE_battery = SE_LiFePO4;
        bat_prefix = strcat('LiFePO4');
        tau_reserve = 0.20;
    case 2
        SE_LiPo = Weight_tier.M_ENERGY.SE_LiPo;
        % Selects the Type of battery
        SE_battery = SE_LiPo;
        bat_prefix = strcat('LiPo');
        tau_reserve = 0.20;
    case 3
        SE_FuelCells = Weight_tier.M_ENERGY.SE_FuelCells;
        % Selects the Type of battery
        SE_battery = SE_FuelCells;
        bat_prefix = strcat('FuelCells');
        tau_reserve = 0.00;
end

% Eficiencias
eta_gear = Prop_data.eta_gear;
eta_m = Prop_data.eta_m;
eta_esc = Prop_data.eta_esc;
eta_dist = Prop_data.eta_dist;
% Initialize Figs
Fig = Fig + 1;

%% Plot options
PERFORMANCE_sensitivity = PLOTS_Prop_Performance.PERFORMANCE_sensitivity;
plots_contour_sensitivity = PLOTS_Prop_Performance.plots_contour_sensitivity;
plots_3d_sensitivity= PLOTS_Prop_Performance.plots_3d_sensitivity;
special_PLOTS = PLOTS_Prop_Performance.special_PLOTS;

for j = 1:length(D_prop_vect)
    D_prop = D_prop_vect(j);
    D_inch = D_prop/2.54*100;
    i = i+1;
    dist = Performance.Range;    %%% Distancia de la misión
 
    %% Velocidad inicial en m/s
    Vh(1) = Performance.V_min;
    
    %% Herramientas bucle
    paro = 0;
    paso = 0.1;
    cont = 1;
    D_inch = D_prop/2.54*100;
    rpm_max = SF_prop*RPMMAX_APC/D_inch;
    nps_max = rpm_max/60;
 
    while paro == 0
        
        % Time of the horizontal flight
        time(cont) = dist/Vh(cont);
        %%% Resistencia aerodinámica
        q_inf = 0.5*rho*Vh(cont)^2;
        CL(cont) = (m_TOW*g)/(q_inf*S_ref);
        CD0 = Aero_TH.CD0;
        CD1 = Aero_TH.CD1;
        CD2 = Aero_TH.CD2;
        CD(cont) = CD0 + CD1*CL(cont) + CD2*CL(cont)^2;
        Dh(cont) = q_inf*S_ref*CD(cont);
        Fdes = Dh(cont); % Matched Drag
        nps_h(cont) = Selection_J_CT_F_design(Vh(cont),Prop_data,Fdes,rho,D_prop,solve_w_fero); % Solves for engine revolucionts for matched drag
        RPM(cont) = nps_h(cont)*60;
        
        % Solves for engine perormance for given Velocity and RPM
        Propulsion_h{cont} = get_EngineProperties_v1(Vh(cont),rho,nps_h(cont),Prop_data,AC_CONFIGURATION,D_prop,n_eng);
        RPM_plot_design_h(cont) = Propulsion_h{cont}.RPM;
            % Coeficientes adimensionales
        CT(cont) = Propulsion_h{cont}.CT;
        CP(cont) = Propulsion_h{cont}.CP;
        CQ(cont) = Propulsion_h{cont}.CQ;
        % eficiencia propulsiva
        eta_p(cont) = Propulsion_h{cont}.etha_ep;
        % eficiencia moto-propulsiva
        eta_total(cont) = Propulsion_h{cont}.etha_emp;
 
        % Total Thrust
        Ti(cont) = Propulsion_h{cont}.Ti;
        % Thrust per engine
        Ti_eng(cont) = Propulsion_h{cont}.Ti_eng;
        % Total Mechanic Power
        Pi(cont) = Propulsion_h{cont}.Pi;
        % Mechanic Power per engine
        Pi_eng(cont) = Propulsion_h{cont}.Pi_eng;
        % Torque per engine
        Qi_eng(cont) = Propulsion_h{cont}.Qi_eng;
        % Total Electri Power
        Pe(cont) = Propulsion_h{cont}.Pe;
        % Electric Power per engine
        Pe_eng(cont) = Propulsion_h{cont}.Pe_eng;
        
        % Advance ratio
        J(cont) = Propulsion_h{cont}.J;
        % Throttle
        delta_Tinter(cont) = nps_h(cont)/nps_max;

        % Energy (W-h)
        Energy(cont) = (1/(1-tau_reserve))*Pe(cont)*time(cont)/3600; % Wh
        Energy_kJ(cont) = (1/(1-tau_reserve))*Pe(cont)*time(cont)/1000; % kJ
        Battery_mass(cont) = Energy(cont)/(SE_battery);
        
        if delta_Tinter(cont)>1
            paro = 1;
            Vmax(j)=Vh(cont);
            nps_graf(j)=nps_h(cont);
        else
            %%% Actualizo variables del bucle
            cont = cont+1;
            Vh(cont) = Vh(cont-1) + paso;
        end
    end
    
    % Thrust
    T_graf(j,1:cont-1) = Ti(1:cont-1);
    T_min(i) = min(Ti);
    T_max(i) = max(Ti);

    % Thrust per engine
    T_graf_eng(j,1:cont-1) = Ti_eng(1:cont-1);
    T_eng_min(i) = min(Ti_eng);
    T_eng_max(i) = max(Ti_eng);

    % Mechanica Power - Climb
    Pi_graf(j,1:cont-1) = Pi(1:cont-1);
    Pi_min(i) = min(Pi);
    Pi_max(i) = max(Pi);

    % Mechanica Power per engine - Climb 
    Pi_eng_graf(j,1:cont-1) = Pi_eng(1:cont-1);
    Pi_eng_min(i) = min(Pi_eng);
    Pi_eng_max(i) = max(Pi_eng);

    % Electric Power - Climb
    Pe_graf(j,1:cont-1) = Pe(1:cont-1);
    Pe_min(i) = min(Pe);
    Pe_max(i) = max(Pe);

    % Electric Power per engine - Climb 
    Pe_eng_graf(j,1:cont-1) = Pe_eng(1:cont-1);
    Pe_eng_min(i) = min(Pe_eng);
    Pe_eng_max(i) = max(Pe_eng);
    
    % Energy (W-h)
    Energy_graf(j,1:cont-1) = Energy(1:cont-1); % Wh
    Energy_min(i) = min(Energy);
    Energy_max(i) = max(Energy);
  
    % Energy (W-h)
    Energy_kJ_graf(j,1:cont-1) = Energy_kJ(1:cont-1); % Wh
    Energy_kJ_min(i) = min(Energy_kJ);
    Energy_kJ_max(i) = max(Energy_kJ);
    
    % Battery mass
    Battery_mass_graf(j,1:cont-1) = Battery_mass(1:cont-1); % kg
    Battery_mass_min(i) = min(Battery_mass);
    Battery_mass_max(i) = max(Battery_mass);
    
    % Torque
    Q_graf(j,1:cont-1) = n_eng*Qi_eng(1:cont-1);
    Q_graf_eng(j,1:cont-1) = Qi_eng(1:cont-1);
    Q_eng_min(i) = min(Qi_eng);
    Q_eng_max(i) = max(Qi_eng);
    
    % Throttle
    delta_T_graf(j,1:cont-1) = delta_Tinter(1:cont-1);
    delta_T_min(i) = min(delta_Tinter);
    delta_T_max(i) = max(delta_Tinter);
    
    % RPMs
    RPM_graf(j,1:cont-1) =  RPM(1:cont-1);
    RPM_min(i) = min(RPM);
    RPM_max(i) = max(RPM);
    
    % Moto-Propulsive Efficiency
    etha_graf(j,1:cont-1) = eta_total(1:cont-1);
    etha_min(i) = min(eta_total);
    etha_max(i) = max(eta_total);
    
    %% Determino el minimo consumo y las nps a las que este se produce
    Th_min(i) = min(Dh);
    [w,z] = find(Dh==Th_min(i));
    nps_dis(i) = nps_h(z);
    Vh_dis(i) = Vh(z);
    [x,y] = find(Energy==Energy_min(i));
    nps_min(i) = nps_h(y);
    Vh_min(i) = Vh(y);
    etha_min(i) = eta_total(y); 
    
end

% Generalized X and Y vector for 3D Plots
X = D_prop_vect*100/2.54; % Prop diameter (inches)
Y = Vh(1:cont-1); % Horizontal Speed
% Slects value for hover as a reference
j_sel = round(length(D_prop_vect)/2);

% Thrust - N
Z_T = T_graf';
vect_cc_T = linspace(min(T_min),max(T_max),N_contour_lines);

% Thrust - N
Z_T_eng = T_graf_eng';
vect_cc_T_eng = linspace(min(T_eng_min),max(T_eng_max),N_contour_lines);

% Mechanica Power - kW
Z_Pi = Pi_graf'/1000;
vect_cc_Pi = linspace(min((Pi_min)/1000),max(Pi_max)/1000,N_contour_lines);

% Mechanica Power per engine - kW
Z_Pi_eng = Pi_eng_graf'/1000;
vect_cc_Pi_eng = linspace(min((Pi_eng_min)/1000),max(Pi_eng_max)/1000,N_contour_lines);

% Electric Power - kW
Z_Pe = Pe_graf'/1000;
vect_cc_Pe = linspace(min((Pe_min)/1000),max(Pe_max)/1000,N_contour_lines);

% Electric Power per engine - kW
Z_Pe_eng = Pe_eng_graf'/1000;
vect_cc_Pe_eng = linspace(min((Pe_eng_min)/1000),max(Pe_eng_max)/1000,N_contour_lines);

% Electric Power per engine - KW-h
Z_Energy = Energy_graf'/1000;
vect_cc_Energy = linspace(min((Energy_min)/1000),max(Energy_max)/1000,N_contour_lines);

% Electric Power per engine  - KW-h
Z_Energy_kJ = Energy_kJ_graf'/1000;
vect_cc_Energy_kJ = linspace(min((Energy_kJ_min)),max(Energy_kJ_max),N_contour_lines);

% delta_T
Z_delta = delta_T_graf';
vect_cc_delta = linspace(min(delta_T_min),max(delta_T_max),N_contour_lines);

% Battery Mass - kg
Z_Battery_mass = Battery_mass_graf';
vect_cc_Battery_mass = linspace(min(Battery_mass_min),max(Battery_mass_max),N_contour_lines);

% RPMs
Z_RPM = RPM_graf';
vect_cc_RPM = linspace(min(RPM_min),max(RPM_max),N_contour_lines);

% Torque
Z_Q = Q_graf_eng';
vect_cc_Q = linspace(min(Q_eng_min),max(Q_eng_max),N_contour_lines);

% Motor-Propulsive Efficiency
Z_etha = etha_graf';
vect_cc_etha = linspace(min(etha_min),max(etha_max),N_contour_lines);

% Storing data
Z_study_h.Z_T = Z_T;
Z_study_h.Z_T_eng = Z_T_eng;
Z_study_h.Z_Pi = Z_Pi;
Z_study_h.Z_Pi_eng = Z_Pi_eng;
Z_study_h.Z_Pe = Z_Pe;
Z_study_h.Z_Pe_eng = Z_Pe_eng;
Z_study_h.Z_Energy = Z_Energy;
Z_study_h.Z_Energy_kJ = Z_Energy_kJ;
Z_study_h.Z_delta = Z_delta;
Z_study_h.Z_Battery_mass = Z_Battery_mass;
Z_study_h.Z_RPM = Z_RPM;
Z_study_h.Z_Q = Z_Q;
Z_study_h.Z_etha = Z_etha;

if PLOTS_SENSITIVITY_HORIZONTAL == 1
    if  special_PLOTS == 1
        figure(Fig)
        [aux H1 H2] = plotyy(D_prop_vect*100/2.54,Energy_min/1000,D_prop_vect*100/2.54,nps_min*60);
        title('Minimum Energy Consumption vs. RPM & o D_p_r_o_p (Horizontal Flight)')
        xlabel('D_p_r_o_p (in)')
        ylabel(aux(1),'kWatt-hour^* - (Horizontal Flight)')
        ylabel(aux(2),'RPMs^*')
        grid on
        name1a   = strcat(prefix,'Min_E_vs_RPM_Dprop',bat_prefix,'_n_eng_',num2str(n_eng));
        if SAVE_FIGS==1
            saveas(gcf,name1a,'fig');
            saveas(gca,name1a,'epsc');
            saveas(gcf,name1a,'pdf');
            saveas(gcf,name1a,'bmp');
            saveas(gcf,name1a,'png');
        end
        Fig = Fig + 1;
        
        figure(Fig)
        [aux H1 H2] = plotyy(D_prop_vect*100/2.54,Th_min,D_prop_vect*100/2.54,nps_min*60);
        title('Minimum Thrust vs. RPM & D_p_r_o_p (Horizontal Flight)')
        xlabel('D_p_r_o_p (in)')
        ylabel(aux(1),'T_h^* (N) - (Horizontal Flight)')
        ylabel(aux(2),'RPMs^*')
        grid on
        name2a   = strcat(prefix,'Min_T_vs_RPM_Dprop',bat_prefix,'_n_eng_',num2str(n_eng));
        if SAVE_FIGS==1
            saveas(gcf,name2a,'fig');
            saveas(gca,name2a,'epsc');
            saveas(gcf,name2a,'pdf');
            saveas(gcf,name2a,'bmp');
            saveas(gcf,name2a,'png');
        end
        Fig = Fig + 1;
        
        figure(Fig)
        [aux H1 H2] = plotyy(D_prop_vect*100/2.54,Vmax,D_prop_vect*100/2.54,nps_graf*60);
        title('V_m_a_x vs. RPM & D_p_r_o_p (Horizontal Flight)')
        xlabel('D_p_r_o_p (in) (in)')
        ylabel(aux(1),'V_{max} (m/s) - (Horizontal Flight)')
        ylabel(aux(2),'RPMs_{max}')
        grid on
        name3a   = strcat(prefix,'Vmax_vs_RPM_Dprop',bat_prefix,'_n_eng_',num2str(n_eng));
        if SAVE_FIGS==1
            saveas(gcf,name3a,'fig');
            saveas(gca,name3a,'epsc');
            saveas(gcf,name3a,'pdf');
            saveas(gcf,name3a,'bmp');
            saveas(gcf,name3a,'png');
        end
        Fig = Fig + 1;
        Warning = 'WARNING!!! Code in PAUSE to check results - PRESS Any Key to Continue';
        disp(Warning)
        pause
    end
    
    if plots_3d_sensitivity ==1
        
        % Legends
        name1   = strcat(prefix,'3D_Thrust_vs_Dprop_Vh_',bat_prefix,'_n_eng_',num2str(n_eng));
        name2   = strcat(prefix,'3D_Pe_eng_vs_Dprop_Vh_',bat_prefix,'_n_eng_',num2str(n_eng));
        name3   = strcat(prefix,'3D_Pe_vs_Dprop_Vh_',bat_prefix,'_n_eng_',num2str(n_eng));
        name4   = strcat(prefix,'3D_deltaT_vs_Dprop_Vh_',bat_prefix,'_n_eng_',num2str(n_eng));
        name5   = strcat(prefix,'3D_RPM_vs_Dprop_Vh_',bat_prefix,'_n_eng_',num2str(n_eng));
        name6   = strcat(prefix,'3D_Q_vs_Dprop_Vh_',bat_prefix,'_n_eng_',num2str(n_eng));
        name7   = strcat(prefix,'3D_E_vs_Dprop_Vh_',bat_prefix,'_n_eng_',num2str(n_eng));
        name7b   = strcat(prefix,'3D_E_CH_vs_Dprop_Vh_',bat_prefix,'_n_eng_',num2str(n_eng));
        name8   = strcat(prefix,'3D_eta_mp_vs_Dprop_Vh_',bat_prefix,'_n_eng_',num2str(n_eng));
        name9   = strcat(prefix,'3D_m_battery_vs_Dprop_Vh_',bat_prefix,'_n_eng_',num2str(n_eng));
        name9b   = strcat(prefix,'3D_m_battery_CH_vs_Dprop_Vh_',bat_prefix,'_n_eng_',num2str(n_eng));
        
        figure(Fig)
        mesh(X,Y,Z_T_eng)
        grid on
        colorbar
        Title1 = strcat('Thrust per Engine vs. D_p_r_o_p & Horizontal Speed');
        title(Title1)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_h (m/s)')
        zlabel('Thrust (N)')
        if SAVE_FIGS==1
            saveas(gcf,name1,'fig');
            saveas(gca,name1,'epsc');
            saveas(gcf,name1,'pdf');
            saveas(gcf,name1,'bmp');
            saveas(gcf,name1,'png');
        end
        Fig = Fig + 1;
        
        figure(Fig)
        mesh(X,Y,Z_Pe_eng)
        grid on
        colorbar
        Title2 = strcat('Electric Power per Engine vs. D_p_r_o_p & Horizontal Speed');
        title(Title2)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_h (m/s)')
        zlabel('P_{eng} (kW)')
        if SAVE_FIGS==1
            saveas(gcf,name2,'fig');
            saveas(gca,name2,'epsc');
            saveas(gcf,name2,'pdf');
            saveas(gcf,name2,'bmp');
            saveas(gcf,name2,'png');
        end
        Fig = Fig + 1;

        figure(Fig)
        mesh(X,Y,Z_Pe)
        grid on
        colorbar
        Title3 = strcat('Total Electric Power vs. D_p_r_o_p & Horizontal Speed');
        title(Title3)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_h (m/s)')
        zlabel('Total P (kW)')
        if SAVE_FIGS==1
            saveas(gcf,name3,'fig');
            saveas(gca,name3,'epsc');
            saveas(gcf,name3,'pdf');
            saveas(gcf,name3,'bmp');
            saveas(gcf,name3,'png');
        end
        Fig = Fig + 1;
        
        figure(Fig)
        mesh(X,Y,Z_delta)
        grid on
        colorbar
        Title4 = strcat('\delta_T per engine vs. D_p_r_o_p & Horizontal Speed');
        title(Title4)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_h (m/s)')
        zlabel('\delta_T')
        if SAVE_FIGS==1
            saveas(gcf,name4,'fig');
            saveas(gca,name4,'epsc');
            saveas(gcf,name4,'pdf');
            saveas(gcf,name4,'bmp');
            saveas(gcf,name4,'png');
        end
        Fig = Fig + 1;
        
        figure(Fig)
        mesh(X,Y,Z_RPM)
        grid on
        colorbar
        Title5 = strcat('RPMs per engine vs. D_p_r_o_p & Horizontal Speed with');
        title(Title5)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_h (m/s)')
        zlabel('RPMs')
        if SAVE_FIGS==1
            saveas(gcf,name5,'fig');
            saveas(gca,name5,'epsc');
            saveas(gcf,name5,'pdf');
            saveas(gcf,name5,'bmp');
            saveas(gcf,name5,'png');
        end
        Fig = Fig + 1;
        
        figure(Fig)
        mesh(X,Y,Z_Q)
        grid on
        colorbar
        Title6 = strcat('Torque (per engine) vs. D_p_r_o_p & Horizontal Speed');
        title(Title6)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_h (m/s)')
        zlabel('Torque Q (Nm)')
        if SAVE_FIGS==1
            saveas(gcf,name6,'fig');
            saveas(gca,name6,'epsc');
            saveas(gcf,name6,'pdf');
            saveas(gcf,name6,'bmp');
            saveas(gcf,name6,'png');
        end
        Fig = Fig + 1;
        
        figure(Fig)
        mesh(X,Y,Z_Energy)
        grid on
        colorbar
        Title7 = strcat('Total Energy vs. D_p_r_o_p & Horizontal Speed');
        title(Title7)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_h (m/s)')
        zlabel('Total Energy (kW-h)')
        if SAVE_FIGS==1
            saveas(gcf,name7,'fig');
            saveas(gca,name7,'epsc');
            saveas(gcf,name7,'pdf');
            saveas(gcf,name7,'bmp');
            saveas(gcf,name7,'png');
        end
        Fig = Fig + 1;

        figure(Fig)
        mesh(X,Y,Z_etha)
        grid on
        colorbar
        Title8 = strcat('Moto-Propulsive Efficiency (\eta_{MP}) vs. D_p_r_o_p & Horizontal Speed');
        title(Title8)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_h (m/s)')
        zlabel('\eta_{MP}')
        if SAVE_FIGS==1
            saveas(gcf,name8,'fig');
            saveas(gca,name8,'epsc');
            saveas(gcf,name8,'pdf');
            saveas(gcf,name8,'bmp');
            saveas(gcf,name8,'png');
        end
        Fig = Fig + 1;

        figure(Fig)
        mesh(X,Y,Z_Battery_mass)
        grid on
        colorbar
        Title9 = strcat('Battery Mass vs. D_p_r_o_p & Horizontal Speed');
        title(Title9)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_h (m/s)')
        zlabel('Battery Mass (kg)')
        if SAVE_FIGS==1
            saveas(gcf,name9,'fig');
            saveas(gca,name9,'epsc');
            saveas(gcf,name9,'pdf');
            saveas(gcf,name9,'bmp');
            saveas(gcf,name9,'png');
        end
        Fig = Fig + 1;
        Warning = 'WARNING!!! Code in PAUSE to check results - PRESS Any Key to Continue';
        disp(Warning)
        pause
    end
    
    % Prints the Contour plots
    if  plots_contour_sensitivity==1
        
        % Legends
        name1   = strcat(prefix,'Thrust_vs_Dprop_Vh_',bat_prefix,'_n_eng_',num2str(n_eng));
        name2   = strcat(prefix,'Pe_eng_vs_Dprop_Vh_',bat_prefix,'_n_eng_',num2str(n_eng));
        name3   = strcat(prefix,'Pe_vs_Dprop_Vh_',bat_prefix,'_n_eng_',num2str(n_eng));
        name4   = strcat(prefix,'deltaT_vs_Dprop_Vh_',bat_prefix,'_n_eng_',num2str(n_eng));
        name5   = strcat(prefix,'RPM_vs_Dprop_Vh_',bat_prefix,'_n_eng_',num2str(n_eng));
        name6   = strcat(prefix,'Q_vs_Dprop_Vh_',bat_prefix,'_n_eng_',num2str(n_eng));
        name7   = strcat(prefix,'E_vs_Dprop_Vh_',bat_prefix,'_n_eng_',num2str(n_eng));
        name7b   = strcat(prefix,'E_CH_vs_Dprop_Vh_',bat_prefix,'_n_eng_',num2str(n_eng));
        name8   = strcat(prefix,'eta_mp_vs_Dprop_Vh_',bat_prefix,'_n_eng_',num2str(n_eng));
        name9   = strcat(prefix,'m_battery_vs_Dprop_Vh_',bat_prefix,'_n_eng_',num2str(n_eng));
        name9b   = strcat(prefix,'m_battery_CH_vs_Dprop_Vh_',bat_prefix,'_n_eng_',num2str(n_eng));
        
        figure(Fig)
        [C_T_c_eng,h_T_c_eng] = contour(X,Y,Z_T_eng,vect_cc_T_eng');
        clabel(C_T_c_eng,h_T_c_eng)
        grid on
        Title1 = strcat('Thrust (N) per Engine vs. D_p_r_o_p & Horizontal Speed');
        title(Title1)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_h (m/s)')
        if SAVE_FIGS==1
            saveas(gcf,name1,'fig');
            saveas(gca,name1,'epsc');
            saveas(gcf,name1,'pdf');
            saveas(gcf,name1,'bmp');
            saveas(gcf,name1,'png');
        end
        
        Fig = Fig + 1;
        
        figure(Fig)
        [C_Pe_eng_c,h_Pe_eng_c] = contour(X,Y,Z_Pe_eng,vect_cc_Pe_eng');
        clabel(C_Pe_eng_c,h_Pe_eng_c)
        grid on
        Title2 = strcat('Electric Power (kW) per Engine vs. D_p_r_o_p & Horizontal Speed');
        title(Title2)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_h (m/s)')
        if SAVE_FIGS==1
            saveas(gcf,name2,'fig');
            saveas(gca,name2,'epsc');
            saveas(gcf,name2,'pdf');
            saveas(gcf,name2,'bmp');
            saveas(gcf,name2,'png');
        end
        Fig = Fig + 1;
        
        figure(Fig)
        [C_Pe_c,h_Pe_c] = contour(X,Y,Z_Pe,vect_cc_Pe');
        clabel(C_Pe_c,h_Pe_c)
        grid on
        Title3 = strcat('Electric Power (kW) vs. D_p_r_o_p & Horizontal Speed');
        title(Title3)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_h (m/s)')
        if SAVE_FIGS==1
            saveas(gcf,name3,'fig');
            saveas(gca,name3,'epsc');
            saveas(gcf,name3,'pdf');
            saveas(gcf,name3,'bmp');
            saveas(gcf,name3,'png');
        end
        Fig = Fig + 1;
        
        figure(Fig)
        [C_delta_c,h_delta_c] = contour(X,Y,Z_delta,vect_cc_delta');
        clabel(C_delta_c,h_delta_c)
        grid on
        Title4 = strcat('\delta_T per engine vs. D_p_r_o_p & Horizontal Speed');
        title(Title4)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_h (m/s)')
        if SAVE_FIGS==1
            saveas(gcf,name4,'fig');
            saveas(gca,name4,'epsc');
            saveas(gcf,name4,'pdf');
            saveas(gcf,name4,'bmp');
            saveas(gcf,name4,'png');
        end
        Fig = Fig + 1;
        
        figure(Fig)
        [C_RPM_c,h_RPM_c] = contour(X,Y,Z_RPM,vect_cc_RPM');
        clabel(C_RPM_c,h_RPM_c)
        grid on
        Title5 = strcat('RPMs per engine vs. D_p_r_o_p & Horizontal Speed');
        title(Title5)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_h (m/s)')
        if SAVE_FIGS==1
            saveas(gcf,name5,'fig');
            saveas(gca,name5,'epsc');
            saveas(gcf,name5,'pdf');
            saveas(gcf,name5,'bmp');
            saveas(gcf,name5,'png');
        end
        Fig = Fig + 1;
        
        figure(Fig)
        [C_Q_c,h_Q_c] = contour(X,Y,Z_Q,vect_cc_Q');
        clabel(C_Q_c,h_Q_c)
        grid on
        Title6 = strcat('Torque (N-m) vs. D_p_r_o_p & Horizontal Speed');
        title(Title6)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_h (m/s)')
        if SAVE_FIGS==1
            saveas(gcf,name6,'fig');
            saveas(gca,name6,'epsc');
            saveas(gcf,name6,'pdf');
            saveas(gcf,name6,'bmp');
            saveas(gcf,name6,'png');
        end
        Fig = Fig + 1;
        
        figure(Fig)
        [C_Energy_c,h_Energy_c] = contour(X,Y,Z_Energy,vect_cc_Energy');
        clabel(C_Energy_c,h_Energy_c)
        grid on
        Title7 = strcat('Total Energy (kW-h) vs. D_p_r_o_p & Horizontal Speed');
        title(Title7)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_h (m/s)')
        if SAVE_FIGS==1
            saveas(gcf,name7,'fig');
            saveas(gca,name7,'epsc');
            saveas(gcf,name7,'pdf');
            saveas(gcf,name7,'bmp');
            saveas(gcf,name7,'png');
        end
        Fig = Fig + 1;
        
        figure(Fig)
        [C_etha_c,h_etha_c] = contour(X,Y,Z_etha,vect_cc_etha');
        clabel(C_etha_c,h_etha_c)
        grid on
        Title8 = strcat('Moto-Propulsive Efficiency (\eta_{MP}) vs. D_p_r_o_p & Horizontal Speed');
        title(Title8)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_h (m/s)')
        if SAVE_FIGS==1
            saveas(gcf,name8,'fig');
            saveas(gca,name8,'epsc');
            saveas(gcf,name8,'pdf');
            saveas(gcf,name8,'bmp');
            saveas(gcf,name8,'png');
        end
        Fig = Fig + 1;
        
        figure(Fig)
        [C_Battery_mass_c,h_Battery_mass_c] = contour(X,Y,Z_Battery_mass,vect_cc_Battery_mass');
        clabel(C_Battery_mass_c,h_Battery_mass_c)
        grid on
        Title9 = strcat('Battery Mass (kg) vs. D_p_r_o_p & Horizontal Speed');
        title(Title9)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_h (m/s)')
        if SAVE_FIGS==1
            saveas(gcf,name9,'fig');
            saveas(gca,name9,'epsc');
            saveas(gcf,name9,'pdf');
            saveas(gcf,name9,'bmp');
            saveas(gcf,name9,'png');
        end
        Fig = Fig + 1;
    end
end