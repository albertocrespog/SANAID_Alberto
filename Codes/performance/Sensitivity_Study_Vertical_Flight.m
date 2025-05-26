function [X,Y,Z_study_v,Fig] = Sensitivity_Study_Vertical_Flight(Aero_TH,SF_prop,RPMMAX_APC,Prop_data,...
    Weight_tier,conv_UNITS,Geo_tier,AC_CONFIGURATION,Fig,PLOTS_SENSITIVITY_VERTICAL,Performance,solve_w_fero,...
    N_prop,pp_D_prop_min,pp_D_prop_max,N_contour_lines,Vv_m,PLOTS_Prop_Performance,Prop_selection,prefix,Plot_Options)

i = 0;
D_prop = Prop_data.D_prop;
S_ref = Geo_tier.S_ref;
rho = Performance.rho;
g = conv_UNITS.g;
m_TOW = Weight_tier.m_TOW;
W_MTOW = m_TOW*g;

SAVE_FIGS = Plot_Options.SAVE_FIGS;
MATLAB_in = Plot_Options.MATLAB_in;
    
n_eng = AC_CONFIGURATION.n_eng;  
D_prop_vect=linspace(pp_D_prop_min*D_prop,pp_D_prop_max*D_prop,N_prop);
dist = Performance.h_climb;
time_Hov = Performance.Endurance_v*60;

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
Hover_PLOTS = PLOTS_Prop_Performance.Hover_PLOTS;
PERFORMANCE_sensitivity = PLOTS_Prop_Performance.PERFORMANCE_sensitivity;
plots_contour_sensitivity = PLOTS_Prop_Performance.plots_contour_sensitivity;
plots_3d_sensitivity= PLOTS_Prop_Performance.plots_3d_sensitivity;
special_PLOTS = PLOTS_Prop_Performance.special_PLOTS;
variation_mass = PLOTS_Prop_Performance.variation_mass;
        
for j = 1:length(D_prop_vect)
    D_prop = D_prop_vect(j);
    D_inch = D_prop/2.54*100;
    i = i+1;
%     dist = Performance.h_climb;
    %% Velocidad inicial en m/s
    Vv(1) = Vv_m;
        
    %%% Herramientas bucle    
    paro = 0;
    paso = 0.1;
    cont = 1;
    D_inch = D_prop/2.54*100;
    rpm_max = SF_prop*RPMMAX_APC/D_inch;
    nps_max = rpm_max/60;

    % Calculates conditions a Hoovering conditions
    if Performance.Endurance_v>0
        V_Hov = 0;
        % Estimation of Desired Thrust
        m_TOW = Weight_tier.m_TOW;
        Fdes_Hov = (m_TOW*g); % For hover condition
        % Function that determines the rev per second for a given prop
        % diameter
        Propulsion_Hov{j} = get_Propulsive_Model_VFP(V_Hov,rho,Fdes_Hov,Prop_data,D_prop_vect(j),n_eng);
%         n_Hov(j) = Selection_J_CT_F_design(V_Hov,Prop_data,Fdes_Hov,rho,D_prop_vect(j),solve_w_fero);
        n_Hov(j) = Propulsion_Hov{j}.RPS;
        delta_T_Hov(j) = n_Hov(j)/nps_max;
        RPM_plot_design_Hov(j) = Propulsion_Hov{j}.RPM;
        T_eng_plot_design_Hov(j) = Propulsion_Hov{j}.Ti_eng;
        Qi_eng_Hov(j) = Propulsion_Hov{j}.Qi_eng;
        Pe_plot_design_Hov(j) = Propulsion_Hov{j}.Pe;
        Pe_eng_plot_design_Hov(j) = Propulsion_Hov{j}.Pe_eng;
        RPM_Hov(j) = n_Hov(j)*60;
%         Propulsion_Hov{j} = get_EngineProperties_v1(V_Hov,rho,n_Hov(j),Prop_data,AC_CONFIGURATION,D_prop_vect(j),n_eng);
%         RPM_plot_design_Hov(j) = Propulsion_Hov{j}.RPM;
%         T_eng_plot_design_Hov(j) = Propulsion_Hov{j}.Ti_eng;
%         Qi_eng_Hov(j) = Propulsion_Hov{j}.Qi_eng;
%         Pe_plot_design_Hov(j) = Propulsion_Hov{j}.Pe;
%         Pe_eng_plot_design_Hov(j) = Propulsion_Hov{j}.Pe_eng;
        Pi_plot_design_Hov(j) = Propulsion_Hov{j}.Pi;
        Pi_eng_plot_design_Hov(j) = Propulsion_Hov{j}.Pi_eng;
%         eta_ES   = eta_m*eta_esc*eta_dist*eta_gear; % Electric system efficiency

%         etha_mp_plot_design_Hov(j) = Propulsion_Hov{j}.ethamp;
%         etha_mp_total_plot_design_Hov(j) = Propulsion_Hov{j}.etha_emp;
        time_Hov(j) = Performance.Endurance_v*60;
        Energy_plot_design_Hov(j) = (1/(1-tau_reserve))*Propulsion_Hov{j}.Pe*time_Hov(j)/3600;  % en W-h
        Energy_plot_design_Hov_kJ(j) = (1/(1-tau_reserve))*Propulsion_Hov{j}.Pe*time_Hov(j)/1000;  % en kJ
        Battery_mass_Hov(j) = Energy_plot_design_Hov(j)/(SE_battery);
    else
        Energy_eng_plot_design_Hov(j) = 0;
        Battery_mass_Hov(j) = 0;
    end
    
    while paro == 0
        time(cont) = dist/Vv(cont);
        %%% Resistencia aerodinámica
        CDv(cont) = 0.5;
        q_inf = 0.5*rho*Vv(cont)^2;
        Dv(cont) = q_inf*S_ref*CDv(cont) + W_MTOW;
        Fdes = Dv(cont);  % Matched Drag       
        nps_v(cont) = Selection_J_CT_F_design(Vv(cont),Prop_data,Fdes,rho,D_prop,solve_w_fero); % Solves for engine revolucionts for matched drag
        RPM(cont) = nps_v(cont)*60;
        % Solves for engine perormance for given Velocity and RPM
        Propulsion_v{cont} = get_EngineProperties_v1(Vv(cont),rho,nps_v(cont),Prop_data,AC_CONFIGURATION,D_prop,n_eng);
        RPM_plot_design_v(cont) = Propulsion_v{cont}.RPM;
        % Coeficientes adimensionales
        CT(cont) = Propulsion_v{cont}.CT;
        CP(cont) = Propulsion_v{cont}.CP;
        CQ(cont) = Propulsion_v{cont}.CQ;
        % eficienci propulsiva
        eta_p(cont) = Propulsion_v{cont}.etha_ep;
        % eficiencia moto-propulsiva
        eta_total(cont) = Propulsion_v{cont}.etha_emp;
                
        % Total Thrust
        Ti(cont) = Propulsion_v{cont}.Ti;
        % Thrust per engine
        Ti_eng(cont) = Propulsion_v{cont}.Ti_eng;
        % Total Mechanic Power
        Pi(cont) = Propulsion_v{cont}.Pi;
        % Mechanic Power per engine
        Pi_eng(cont) = Propulsion_v{cont}.Pi_eng;
        % Torque per engine
        Qi_eng(cont) = Propulsion_v{cont}.Qi_eng;
        % Total Electri Power
        Pe(cont) = Propulsion_v{cont}.Pe;
        % Electric Power per engine
        Pe_eng(cont) = Propulsion_v{cont}.Pe_eng;

        % Advance ratio
        J(cont) = Propulsion_v{cont}.J;
        % Throttle
        delta_Tinter(cont) = nps_v(cont)/nps_max;

        % Energy (W-h)
        Energy(cont) = (1/(1-tau_reserve))*Pe(cont)*time(cont)/(3600); % Wh
        Energy_kJ(cont) = (1/(1-tau_reserve))*Pe(cont)*time(cont)/1000; % kJ
        Battery_mass(cont) = Energy(cont)/(SE_battery);
        Energy_plot_CH(cont) = Energy(cont) + Energy_plot_design_Hov(j);
        % Combined Energy - climb and hover - kJ
        Energy_plot_CH_kJ(cont) = Energy_kJ(cont) + Energy_plot_design_Hov_kJ(j);
        % Combined battery mass - climb and hover
        Battery_mass_CH(cont) = Battery_mass(cont) + Battery_mass_Hov(j);
        
        if delta_Tinter(cont)>1
            paro = 1;
            Vmax(j)=Vv(cont);
            nps_graf(j)=nps_v(cont);
        else
            % Actualizo variables del bucle
            cont = cont+1;
            Vv(cont) = Vv(cont-1)+paso;
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
    
    % Combined Energy - Climb and hover - Wh
    Energy_CH_graf(j,1:cont-1) = Energy_plot_CH(1:cont-1);
    Energy_CH_min(i) = min(Energy_plot_CH);
    Energy_CH_max(i) = max(Energy_plot_CH);
    
    % Combined Energy - Climb and hover - kJ
    Energy_CH_kJ_graf(j,1:cont-1) = Energy_plot_CH_kJ(1:cont-1);
    Energy_CH_kJ_min(i) = min(Energy_plot_CH_kJ);
    Energy_CH_kJ_max(i) = max(Energy_plot_CH_kJ);
    
    % Combined battery mass - climb and hover
    Battery_mass_CH_graf(j,1:cont-1) = Battery_mass_CH(1:cont-1);
    Battery_mass_CH_min(i) = min(Battery_mass_CH);
    Battery_mass_CH_max(i) = max(Battery_mass_CH);
    
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
    Tv_min(i) = min(Dv);
    [w,z] = find(Dv==Tv_min(i));
    nps_dis(i) = nps_v(z);
    Vv_dis(i) = Vv(z);
    [x,y] = find(Energy==Energy_min(i));
    nps_min(i) = nps_v(y);
    Vv_min(i) = Vv(y);
    etha_min(i) = eta_p(y); 
end

if variation_mass == 1
    [DATA_mass,Fig] = variation_mass_function_v1(Prop_data,D_prop_vect,AC_CONFIGURATION,SF_prop,...
        SE_battery,RPMMAX_APC,Weight_tier,Fig,Prop_selection,N_contour_lines,prefix,bat_prefix,Plot_Options,tau_reserve);
end

% Generalized X and Y vector for 3D Plots
X = D_prop_vect*100/2.54; % Prop diameter (inches)
Y = Vv(1:cont-1); % Horizontal Speed
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

% Combined (climb and hover)Electric Power per engine - KW-h
Z_Energy_CH = Energy_CH_graf'/1000;
vect_cc_Energy_CH = linspace(min((Energy_CH_min)/1000),max(Energy_CH_max)/1000,N_contour_lines);

% Combined (climb and hover)Electric Power per engine  - KW-h
Z_Energy_CH_kJ = Energy_CH_kJ_graf'/1000;
vect_cc_Energy_CH_kJ = linspace(min((Energy_CH_kJ_min)),max(Energy_CH_kJ_max),N_contour_lines);

% Combined (climb and hover) Battery Mass - kg
Z_Battery_mass_CH = Battery_mass_CH_graf';
vect_cc_Battery_mass_CH = linspace(min(Battery_mass_CH_min),max(Battery_mass_CH_max),N_contour_lines);

% Storing data
Z_study_v.Z_T = Z_T;
Z_study_v.Z_T_eng = Z_T_eng;
Z_study_v.Z_Pi = Z_Pi;
Z_study_v.Z_Pi_eng = Z_Pi_eng;
Z_study_v.Z_Pe = Z_Pe;
Z_study_v.Z_Pe_eng = Z_Pe_eng;
Z_study_v.Z_Energy = Z_Energy;
Z_study_v.Z_Energy_kJ = Z_Energy_kJ;
Z_study_v.Z_delta = Z_delta;
Z_study_v.Z_Battery_mass = Z_Battery_mass;
Z_study_v.Z_RPM = Z_RPM;
Z_study_v.Z_Q = Z_Q;
Z_study_v.Z_etha = Z_etha;
Z_study_v.Z_Energy_CH = Z_Energy_CH;
Z_study_v.Z_Energy_CH_kJ = Z_Energy_CH_kJ;
Z_study_v.Z_Battery_mass_CH = Z_Battery_mass_CH;

if PLOTS_SENSITIVITY_VERTICAL == 1
    % Conducts the plots for the hover conditions
    if Hover_PLOTS == 1
    
        % Legend
        name1   = strcat(prefix,'Thrust_vs_Dprop_Hover_',bat_prefix,'_n_eng_',num2str(n_eng));
        name2   = strcat(prefix,'Pe_eng_vs_Dprop_Hover_',bat_prefix,'_n_eng_',num2str(n_eng));
        name3   = strcat(prefix,'Pe_vs_Dprop_Hover_',bat_prefix,'_n_eng_',num2str(n_eng));
        name4   = strcat(prefix,'deltaT_vs_Dprop_Hover_',bat_prefix,'_n_eng_',num2str(n_eng));
        name5   = strcat(prefix,'RPM_vs_Dprop_Hover_',bat_prefix,'_n_eng_',num2str(n_eng));
        name6   = strcat(prefix,'Q_vs_Dprop_Hover_',bat_prefix,'_n_eng_',num2str(n_eng));
        name7   = strcat(prefix,'E_vs_Dprop_Hover_',bat_prefix,'_n_eng_',num2str(n_eng));
        name7b   = strcat(prefix,'E_CH_vs_Dprop_Hover_',bat_prefix,'_n_eng_',num2str(n_eng));
        name8   = strcat(prefix,'eta_mp_vs_Dprop_Hover_',bat_prefix,'_n_eng_',num2str(n_eng));
        name9   = strcat(prefix,'m_battery_vs_Dprop_Hover_',bat_prefix,'_n_eng_',num2str(n_eng));
        name9b   = strcat(prefix,'m_battery_CH_vs_Dprop_Hover_',bat_prefix,'_n_eng_',num2str(n_eng));

        figure(Fig)
        plot(D_prop_vect*100/2.54,T_eng_plot_design_Hov);
        title('Thrust vs. D_p_r_o_p (Hover)')
        xlabel('D_p_r_o_p (in)')
        ylabel('T (N)')
        grid on
        if SAVE_FIGS==1
            saveas(gcf,name1,'fig');
            saveas(gca,name1,'epsc');
            saveas(gcf,name1,'pdf');
            saveas(gcf,name1,'bmp');
            saveas(gcf,name1,'png');
        end
        Fig = Fig + 1;

        figure(Fig)
        plot(D_prop_vect*100/2.54,Qi_eng_Hov);
        title('Torque vs. D_p_r_o_p (Hover)')
        xlabel('D_p_r_o_p (in)')
        ylabel('Q (N-m)')
        grid on
        if SAVE_FIGS==1
            saveas(gcf,name6,'fig');
            saveas(gca,name6,'epsc');
            saveas(gcf,name6,'pdf');
            saveas(gcf,name6,'bmp');
            saveas(gcf,name6,'png');
        end
        Fig = Fig + 1;

        figure(Fig)
        plot(D_prop_vect*100/2.54,Pe_eng_plot_design_Hov/1000);
        title('Electric Power per engine vs. D_p_r_o_p (Hover)')
        xlabel('D_p_r_o_p (in)')
        ylabel('P_e per engine(kW)')
        grid on
        if SAVE_FIGS==1
            saveas(gcf,name2,'fig');
            saveas(gca,name2,'epsc');
            saveas(gcf,name2,'pdf');
            saveas(gcf,name2,'bmp');
            saveas(gcf,name2,'png');
        end
        Fig = Fig + 1;

        figure(Fig)
        plot(D_prop_vect*100/2.54,Energy_plot_design_Hov/1000);
        title('Energy vs. D_p_r_o_p (Hover)')
        xlabel('D_p_r_o_p (in)')
        ylabel('Energy(kWh)')
        grid on
        if SAVE_FIGS==1
            saveas(gcf,name7,'fig');
            saveas(gca,name7,'epsc');
            saveas(gcf,name7,'pdf');
            saveas(gcf,name7,'bmp');
            saveas(gcf,name7,'png');
        end

        Fig = Fig + 1;

        figure(Fig)
        plot(D_prop_vect*100/2.54,Battery_mass_Hov);
        title('Battery Mass per engine vs. D_p_r_o_p (Hover)')
        xlabel('D_p_r_o_p (in)')
        ylabel('m_{bat} (kg)')
        grid on
        Fig = Fig + 1;
        if SAVE_FIGS==1
            saveas(gcf,name9,'fig');
            saveas(gca,name9,'epsc');
            saveas(gcf,name9,'pdf');
            saveas(gcf,name9,'bmp');
            saveas(gcf,name9,'png');
        end
        Warning = 'WARNING!!! Code in PAUSE to check results - PRESS Any Key to Continue';
        disp(Warning)
        pause
    end
    
    if  special_PLOTS == 1
        figure(Fig)
        [aux H1 H2] = plotyy(D_prop_vect*100/2.54,Energy_min/1000,D_prop_vect*100/2.54,nps_min*60);
        title('Minimum Energy Consumption vs. RPM & o D_p_r_o_p (Axial Flight)')
        xlabel('D_p_r_o_p')
        ylabel(aux(1),'Energy^* (Axial Flight) (kW-h)')
        ylabel(aux(2),'RPMs^* (Axial Flight)')
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
        [aux H1 H2] = plotyy(D_prop_vect*100/2.54,Tv_min,D_prop_vect*100/2.54,nps_min*60);
        title('Minimum Thrust vs. RPM & D_p_r_o_p (Axial Flight)')
        xlabel('D_p_r_o_p')
        ylabel(aux(1),'T_v^* (N)')
        ylabel(aux(2),'RPMs^* (Axial Flight)')
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
        title('V_m_a_x vs. RPM & D_p_r_o_p (Axial Flight)')
        xlabel('D_p_r_o_p (in)')
        ylabel(aux(1),'V_m_a_x (Axial Flight)')
        ylabel(aux(2),'RPMs Max')
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

% Prints the 3D plots
    if plots_3d_sensitivity ==1
        % Legends
        name1   = strcat(prefix,'3D_Thrust_vs_Dprop_Vv_',bat_prefix,'_n_eng_',num2str(n_eng));
        name2   = strcat(prefix,'3D_Pe_eng_vs_Dprop_Vv_',bat_prefix,'_n_eng_',num2str(n_eng));
        name3   = strcat(prefix,'3D_Pe_vs_Dprop_Vv_',bat_prefix,'_n_eng_',num2str(n_eng));
        name4   = strcat(prefix,'3D_deltaT_vs_Dprop_Vv_',bat_prefix,'_n_eng_',num2str(n_eng));
        name5   = strcat(prefix,'3D_RPM_vs_Dprop_Vv_',bat_prefix,'_n_eng_',num2str(n_eng));
        name6   = strcat(prefix,'3D_Q_vs_Dprop_Vv_',bat_prefix,'_n_eng_',num2str(n_eng));
        name7   = strcat(prefix,'3D_E_vs_Dprop_Vv_',bat_prefix,'_n_eng_',num2str(n_eng));
        name7b   = strcat(prefix,'3D_E_CH_vs_Dprop_Vv_',bat_prefix,'_n_eng_',num2str(n_eng));
        name8   = strcat(prefix,'3D_eta_mp_vs_Dprop_Vv_',bat_prefix,'_n_eng_',num2str(n_eng));
        name9   = strcat(prefix,'3D_m_battery_vs_Dprop_Vv_',bat_prefix,'_n_eng_',num2str(n_eng));
        name9b   = strcat(prefix,'3D_m_battery_CH_vs_Dprop_Vv_',bat_prefix,'_n_eng_',num2str(n_eng));
        
        figure(Fig)
        mesh(X,Y,Z_T_eng)
        grid on
        colorbar
        Title1 = strcat('Thrust per Engine vs. D_p_r_o_p & Vertical Speed with T_{Hover}=',num2str(T_eng_plot_design_Hov(j_sel)),' N, @ D_{prop}=',num2str(X(j_sel)),'in');
        title(Title1)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_v (m/s)')
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
        Title2 = strcat('Electric Power per Engine vs. D_p_r_o_p & Vertical Speed with P_{Hover}=',num2str(Pe_eng_plot_design_Hov(j_sel)/1000),' kW, @ D_{prop}=',num2str(X(j_sel)),'in');
        title(Title2)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_v (m/s)')
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
        Title3 = strcat('Total Electric Power vs. D_p_r_o_p & Vertical Speed with P_{Hover}=',num2str(Pe_plot_design_Hov(j_sel)/1000),' kW, @ D_{prop}=',num2str(X(j_sel)),'in');
        title(Title3)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_v (m/s)')
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
        Title4 = strcat('\delta_T per engine vs. D_p_r_o_p & Vertical Speed with \delta_{Hover}=',num2str(delta_T_Hov(j_sel)),' %, @ D_{prop}=',num2str(X(j_sel)),'in');
        title(Title4)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_v (m/s) - (Vertical Flight)')
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
        Title5 = strcat('RPMs per engine vs. D_p_r_o_p & Vertical Speed with RPMs_{Hover}=',num2str(RPM_Hov(j_sel)),', @ D_{prop}=',num2str(X(j_sel)),'in');
        title(Title5)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_v (m/s) - (Vertical Flight)')
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
        Title6 = strcat('Torque (per engine) vs. D_p_r_o_p & Vertical Speed with Q_{Hover}=',num2str(Qi_eng_Hov(j_sel)),' Nm, @ D_{prop}=',num2str(X(j_sel)),'in');
        title(Title6)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_v (m/s) - (Vertical Flight)')
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
        Title7 = strcat('Total Energy vs. D_p_r_o_p & Vertical Speed with E_{Hoovering}=',num2str(Energy_plot_design_Hov(j_sel)/1000),' kW-h, @ D_{prop}=',num2str(X(j_sel)),'in');
        title(Title7)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_v (m/s)')
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
        mesh(X,Y,Z_Energy_CH)
        grid on
        colorbar
        Title7b = strcat('Total Energy (Climb and Hover) vs. D_p_r_o_p & Vertical Speed with E_{Hoovering}=',num2str(Energy_plot_design_Hov(j_sel)/1000),' kW-h, @ D_{prop}=',num2str(X(j_sel)),'in');
        title(Title7b)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_v (m/s)')
        zlabel('Total Energy (kW-h) - (Climb and Hover)')
        if SAVE_FIGS==1
            saveas(gcf,name7b,'fig');
            saveas(gca,name7b,'epsc');
            saveas(gcf,name7b,'pdf');
            saveas(gcf,name7b,'bmp');
            saveas(gcf,name7b,'png');
        end
        Fig = Fig + 1;

        figure(Fig)
        mesh(X,Y,Z_etha)
        grid on
        colorbar
        Title8 = strcat('Moto-Propulsive Efficiency (\eta_{MP}) vs. D_p_r_o_p & Vertical Speed with');
        title(Title8)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_v (m/s)')
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
        Title9 = strcat('Battery Mass vs. D_p_r_o_p & Vertical Speed with m_{bat}=',num2str(Battery_mass_Hov(j_sel)),', @ D_{prop}=',num2str(X(j_sel)),'in');
        title(Title9)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_v (m/s)')
        zlabel('Battery Mass (kg)')
        if SAVE_FIGS==1
            saveas(gcf,name9,'fig');
            saveas(gca,name9,'epsc');
            saveas(gcf,name9,'pdf');
            saveas(gcf,name9,'bmp');
            saveas(gcf,name9,'png');
        end
        Fig = Fig + 1;

        figure(Fig)
        mesh(X,Y,Z_Battery_mass_CH)
        grid on
        colorbar
        Title9 = strcat('Battery Mass (Climb & Hover) vs. D_p_r_o_p & Vertical Speed with m_{bat}=',num2str(Battery_mass_Hov(j_sel)),', @ D_{prop}=',num2str(X(j_sel)),'in');
        title(Title9)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_v (m/s)')
        zlabel('Battery Mass (kg) - (Climb & Hover)')
         if SAVE_FIGS==1
            saveas(gcf,name9b,'fig');
            saveas(gca,name9b,'epsc');
            saveas(gcf,name9b,'pdf');
            saveas(gcf,name9b,'bmp');
            saveas(gcf,name9b,'png');
        end
       Fig = Fig + 1;
       Warning = 'WARNING!!! Code in PAUSE to check results - PRESS Any Key to Continue';
       disp(Warning)
       pause
    end
    
    % Prints the Contour plots
    if  plots_contour_sensitivity==1
        
        % Legends
        name1   = strcat(prefix,'Thrust_vs_Dprop_Vv_',bat_prefix,'_n_eng_',num2str(n_eng));
        name2   = strcat(prefix,'Pe_eng_vs_Dprop_Vv_',bat_prefix,'_n_eng_',num2str(n_eng));
        name3   = strcat(prefix,'Pe_vs_Dprop_Vv_',bat_prefix,'_n_eng_',num2str(n_eng));
        name4   = strcat(prefix,'deltaT_vs_Dprop_Vv_',bat_prefix,'_n_eng_',num2str(n_eng));
        name5   = strcat(prefix,'RPM_vs_Dprop_Vv_',bat_prefix,'_n_eng_',num2str(n_eng));
        name6   = strcat(prefix,'Q_vs_Dprop_Vv_',bat_prefix,'_n_eng_',num2str(n_eng));
        name7   = strcat(prefix,'E_vs_Dprop_Vv_',bat_prefix,'_n_eng_',num2str(n_eng));
        name7b   = strcat(prefix,'E_CH_vs_Dprop_Vv_',bat_prefix,'_n_eng_',num2str(n_eng));
        name8   = strcat(prefix,'eta_mp_vs_Dprop_Vv_',bat_prefix,'_n_eng_',num2str(n_eng));
        name9   = strcat(prefix,'m_battery_vs_Dprop_Vv_',bat_prefix,'_n_eng_',num2str(n_eng));
        name9b   = strcat(prefix,'m_battery_CH_vs_Dprop_Vv_',bat_prefix,'_n_eng_',num2str(n_eng));
        
        figure(Fig)
        [C_T_c_eng,h_T_c_eng] = contour(X,Y,Z_T_eng,vect_cc_T_eng');
        clabel(C_T_c_eng,h_T_c_eng)
        grid on
        Title1 = strcat('Thrust (N) per Engine vs. D_p_r_o_p & Vertical Speed with T_{Hover}=',num2str(T_eng_plot_design_Hov(j_sel)),' N, @ D_{prop}=',num2str(X(j_sel)),'in');
        title(Title1)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_v (m/s)')
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
        Title2 = strcat('Electric Power (kW) per Engine vs. D_p_r_o_p & Vertical Speed with P_{Hover}=',num2str(Pe_eng_plot_design_Hov(j_sel)/1000),' kW, @ D_{prop}=',num2str(X(j_sel)),'in');
        title(Title2)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_v (m/s)')
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
        Title3 = strcat('Electric Power (kW) vs. D_p_r_o_p & Vertical Speed with P_{Hover}=',num2str(Pe_plot_design_Hov(j_sel)/1000),' kW, @ D_{prop}=',num2str(X(j_sel)),'in');
        title(Title3)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_v (m/s)')
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
        Title4 = strcat('\delta_T per engine vs. D_p_r_o_p & Vertical Speed with \delta_{Hover}=',num2str(delta_T_Hov(j_sel)),' %, @ D_{prop}=',num2str(X(j_sel)),'in');
        title(Title4)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_v (m/s)')
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
        Title5 = strcat('RPMs per engine vs. D_p_r_o_p & Vertical Speed with RPMs_{Hover}=',num2str(RPM_Hov(j_sel)),', @ D_{prop}=',num2str(X(j_sel)),'in');
        title(Title5)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_v (m/s)')
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
        Title6 = strcat('Torque (N-m) vs. D_p_r_o_p & Vertical Speed with Q_{Hover}=',num2str(Qi_eng_Hov(j_sel)),' Nm, @ D_{prop}=',num2str(X(j_sel)),'in');
        title(Title6)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_v (m/s)')
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
        Title7 = strcat('Total Energy (kW-h) vs. D_p_r_o_p & Vertical Speed with E_{Hoovering}=',num2str(Energy_plot_design_Hov(j_sel)/1000),' kW-h, @ D_{prop}=',num2str(X(j_sel)),'in');
        title(Title7)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_v (m/s)')
         if SAVE_FIGS==1
            saveas(gcf,name7,'fig');
            saveas(gca,name7,'epsc');
            saveas(gcf,name7,'pdf');
            saveas(gcf,name7,'bmp');
            saveas(gcf,name7,'png');
        end
        Fig = Fig + 1;

        figure(Fig)
        [C_Energy_CH_c,h_Energy_CH_c] = contour(X,Y,Z_Energy_CH,vect_cc_Energy_CH');
        clabel(C_Energy_CH_c,h_Energy_CH_c)
        grid on
        Title7b = strcat('Total Energy (Climb and Hover) (kW-h) vs. D_p_r_o_p & Vertical Speed with E_{Hoovering}=',num2str(Energy_plot_design_Hov(j_sel)/1000),' kW-h, @ D_{prop}=',num2str(X(j_sel)),'in');
        title(Title7b)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_v (m/s)')
         if SAVE_FIGS==1
            saveas(gcf,name7b,'fig');
            saveas(gca,name7b,'epsc');
            saveas(gcf,name7b,'pdf');
            saveas(gcf,name7b,'bmp');
            saveas(gcf,name7b,'png');
        end
        Fig = Fig + 1;

        figure(Fig)
        [C_etha_c,h_etha_c] = contour(X,Y,Z_etha,vect_cc_etha');
        clabel(C_etha_c,h_etha_c)
        grid on
        Title8 = strcat('Moto-Propulsive Efficiency (\eta_{MP}) vs. D_p_r_o_p & Vertical Speed');
        title(Title8)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_v (m/s)')
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
        Title9 = strcat('Battery Mass vs. D_p_r_o_p (kg) & Vertical Speed with m_{bat}=',num2str(Battery_mass_Hov(j_sel)),', @ D_{prop}=',num2str(X(j_sel)),'in');
        title(Title9)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_v (m/s)')
         if SAVE_FIGS==1
            saveas(gcf,name9,'fig');
            saveas(gca,name9,'epsc');
            saveas(gcf,name9,'pdf');
            saveas(gcf,name9,'bmp');
            saveas(gcf,name9,'png');
        end
        Fig = Fig + 1;

        figure(Fig)
        [C_Battery_mass_CH_c,h_Battery_mass_CH_c] = contour(X,Y,Z_Battery_mass_CH,vect_cc_Battery_mass_CH');
        clabel(C_Battery_mass_CH_c,h_Battery_mass_CH_c)
        grid on
        Title9 = strcat('Battery Mass (Climb & Hover) (kg) vs. D_p_r_o_p & Vertical Speed with m_{bat}=',num2str(Battery_mass_Hov(j_sel)),', @ D_{prop}=',num2str(X(j_sel)),'in');
        title(Title9)
        xlabel('D_p_r_o_p (in)')
        ylabel('V_v (m/s)')
         if SAVE_FIGS==1
            saveas(gcf,name9b,'fig');
            saveas(gca,name9b,'epsc');
            saveas(gcf,name9b,'pdf');
            saveas(gcf,name9b,'bmp');
            saveas(gcf,name9b,'png');
        end
        Fig = Fig + 1;
        Warning = 'WARNING!!! Code in PAUSE to check results - PRESS Any Key to Continue';
        disp(Warning)
        pause
    end
end