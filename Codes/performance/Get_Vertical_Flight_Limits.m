function [Fig] = Get_Vertical_Flight_Limits(Aero_TH,pp,RPMMAX_APC,Prop_data,Weight_tier,conv_UNITS,Performance_preliminar,Geo_tier,...
    AC_CONFIGURATION,N_Vv,Vv_vect,Fig,PLOT_Get_Vertical_Flight_Limits,solve_w_fero,Plot_Options)

rho = Performance_preliminar.rho;
dist = Performance_preliminar.h_climb;
MATLAB_in = Plot_Options.MATLAB_in;

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

n_eng = AC_CONFIGURATION.n_eng;  

S_ref = Geo_tier.S_ref;
g = conv_UNITS.g;

% Initialize Figs
Fig = Fig + 1;
% while convergence == 1
Delta_Vv = 0.1;
for i = 1:N_Vv
            Vv = Vv_vect(i);
            D_prop(i) = determination_D_prop_vertical(Aero_TH,pp,RPMMAX_APC,Vv,Prop_data,Weight_tier,...
                conv_UNITS,Performance_preliminar,Geo_tier,AC_CONFIGURATION);
            D_prop_inch(i) = D_prop(i)*100/2.54;
            
            % Estimation of Desired Thrust
            q_inf = 0.5*rho*Vv^2;
            m_TOW = Weight_tier.m_TOW;
            
            C_DV=0.5;
            Dv = 0.5*rho*(Vv^2)*S_ref*C_DV + (m_TOW*g);

            Fdes = Dv;
            n(i) = Selection_J_CT_F_design(Vv,Prop_data,Fdes,rho,D_prop(i),solve_w_fero);
            RPM(i) = n(i)*60;
            
%             Propulsion{i} = get_EngineProperties_design(Vv,rho,n(i),Prop_data,AC_CONFIGURATION,D_prop(i),n_eng);
            Propulsion{i} = get_EngineProperties_v1(Vv,rho,n(i),Prop_data,AC_CONFIGURATION,D_prop(i),n_eng);
             
            RPM_plot_design(i) = Propulsion{i}.RPM;
            T_eng_plot_design(i) = Propulsion{i}.Ti_eng;
            Pe_eng_plot_design(i) = Propulsion{i}.Pe_eng;
            Pe_plot_design(i) = Propulsion{i}.Pe;
            etha_ep_plot_design(i) = Propulsion{i}.etha_ep;
            etha_em_plot_design(i) = Propulsion{i}.etha_em;
            etha_emp_plot_design(i) = Propulsion{i}.etha_emp;
            time_V(i) = dist/Vv;
            Energy_eng_plot_design(i) = (1/(1-tau_reserve))*Propulsion{i}.Pe*time_V(i)/3600; % En W-h 
            Energy_eng_plot_design_kJ(i) = (1/(1-tau_reserve))*Propulsion{i}.Pe*time_V(i)/1000; % En kJ 
%             Battery_mass(i) = Energy_eng_plot_design(i)/(SE_LiFePO4);
            Battery_mass(i) = Energy_eng_plot_design(i)/(SE_battery);
end

% Calculates conditions a Hoovering conditions
if Performance_preliminar.Endurance_v>0 
    
    V_Hov = 0;
    D_prop_Hov = determination_D_prop_vertical(Aero_TH,pp,RPMMAX_APC,V_Hov,Prop_data,Weight_tier,...
        conv_UNITS,Performance_preliminar,Geo_tier,AC_CONFIGURATION);
    D_prop_inch_Hov = D_prop_Hov*100/2.54;
    
%     V_batt =  (pp*RPMMAX_APC/(100*D_prop_Hov/2.54))/195;
%     v_cell = 3.3;
%     num_cells = round(V_batt/v_cell) + 1
%     Weight_batt = num_cells*(76/1000)
%     Capacity = 2300; % 2300 mAh
%     Wbat = (Capacity/1000)*V_batt;
%     P = V_batt*Ic
%     Wb*Ib/Wb*Ic = t
    % Estimation of Desired Thrust
    m_TOW = Weight_tier.m_TOW;
      
    Fdes_Hov = (m_TOW*g);
    n_Hov = Selection_J_CT_F_design(V_Hov,Prop_data,Fdes_Hov,rho,D_prop_Hov,solve_w_fero);
    RPM_Hov = n_Hov*60;
%     Propulsion_Hov = get_EngineProperties_design(V_Hov,rho,n_Hov,Prop_data,AC_CONFIGURATION,D_prop_Hov);
    Propulsion_Hov = get_EngineProperties_v1(V_Hov,rho,n_Hov,Prop_data,AC_CONFIGURATION,D_prop_Hov,n_eng);
    
    RPM_plot_design_Hov = Propulsion_Hov.RPM;
    T_eng_plot_design_Hov = Propulsion_Hov.Ti_eng;
    Pe_eng_plot_design_Hov = Propulsion_Hov.Pe_eng;
    Pe_plot_design_Hov = Propulsion_Hov.Pe;
    etha_em_plot_design_Hov = Propulsion_Hov.etha_em;
    etha_ep_plot_design_Hov = Propulsion_Hov.etha_ep;
    etha_emp_plot_design_Hov = Propulsion_Hov.etha_emp;
    time_Hov = Performance_preliminar.Endurance_v*60;
    Energy_eng_plot_design_Hov = (1/(1-tau_reserve))*Propulsion_Hov.Pe*time_Hov/3600;  % en W-h
    Energy_eng_plot_design_Hov_kJ = (1/(1-tau_reserve))*Propulsion_Hov.Pe*time_Hov/1000;  % en kJ
%     Battery_mass_Hov = Energy_eng_plot_design_Hov/(SE_LiFePO4);  
    Battery_mass_Hov = Energy_eng_plot_design_Hov/(SE_battery);  
else
    Energy_eng_plot_design_Hov = 0;
    Battery_mass_Hov = 0;
end

LS = Plot_Options.LS;
FS = Plot_Options.FS;
LFS = Plot_Options.LFS*0.75;
leg = strcat('W_{payload}= ',num2str(Weight_tier.M_PAYLOAD.m_cargo),'kg & Safety Factor =',num2str(pp));


if PLOT_Get_Vertical_Flight_Limits == 1
    figure(Fig)
    subplot(2,2,1)
    plot(Vv_vect,D_prop_inch)
    title('Prop Diameter vs. Vertical Speed')
    xlabel('Vertical Speed - V_v (m/s)')
    ylabel('Prop Diameter - D (inch)')
    if MATLAB_in == 1
        legend(leg)
    else
        legend(leg,0,'FontSize',LFS)
    end
    grid on
    
    subplot(2,2,2)
    plot(Vv_vect,T_eng_plot_design)
    title('Thrust per Engine vs. Vertical Speed')
    xlabel('Vertical Speed - V_v (m/s)')
    ylabel('Thrust per Engine - (T_{eng}) (N)')
    if MATLAB_in == 1
        legend(leg)
    else
        h_legend = legend(leg)
        set(h_legend, 'Location','Best','FontSize',LFS)
    end
    grid on
    
    subplot(2,2,3)
    plot(Vv_vect,Pe_plot_design/1000)
    Title1 = strcat('Total Electric Power vs. Vertical Speed with P_{hover}=',num2str(Pe_plot_design_Hov/1000),' kW');
    title(Title1)
    xlabel('Vertical Speed - V_v (m/s)')
    ylabel('Total Electric Power - (P)(KW)')
    if MATLAB_in == 1
        legend(leg)
    else
        h_legend = legend(leg)
        set(h_legend, 'Location','Best','FontSize',LFS)
    end
    grid on

    subplot(2,2,4)
    plot(Vv_vect,Pe_eng_plot_design/1000)
    Title2 = strcat('Electric Power per Engine vs. Vertical Speed with P_{hover-eng}=',num2str(Pe_eng_plot_design_Hov/1000),' kW');
    title(Title2)
    xlabel('Vertical Speed - V_v (m/s)')
    ylabel('Electric Power per Engine - (P_{hover-eng})(KW)')
    if MATLAB_in == 1
        legend(leg)
    else
        h_legend = legend(leg)
        set(h_legend, 'Location','Best','FontSize',LFS)
    end
    grid on

    Fig = Fig + 1;
    figure(Fig)
    subplot(2,2,1)
    plot(Vv_vect,(Energy_eng_plot_design)/1000)
    Title3 = strcat('Climb Energy vs. Vertical Speed with E_{hover}=',num2str(Energy_eng_plot_design_Hov/1000),' kWh for t_{hover}=',num2str(time_Hov),'seg');
    title(Title3)
    xlabel('Vertical Speed - V_v (m/s)')
    ylabel('Climb and Hover Energy - E_{climb & hover} (kW-h)')
    if MATLAB_in == 1
        legend(leg)
    else
        h_legend = legend(leg)
        set(h_legend, 'Location','Best','FontSize',LFS)
    end
    grid on

    subplot(2,2,2)
    plot(Vv_vect,(Energy_eng_plot_design+Energy_eng_plot_design_Hov)/1000)
    Title3 = strcat('Climb and Hover Energy vs. Vertical Speed with E_{hover}=',num2str(Energy_eng_plot_design_Hov/1000),' kWh for t_{hover}=',num2str(time_Hov),'seg');
    title(Title3)
    xlabel('Vertical Speed - V_v (m/s)')
    ylabel('Climb and Hover Energy - E_{climb & hover} (kW-h)')
    if MATLAB_in == 1
        legend(leg)
    else
        h_legend = legend(leg)
        set(h_legend, 'Location','Best','FontSize',LFS)
    end
    grid on
    
    subplot(2,2,3)
    plot(Vv_vect,etha_emp_plot_design)
    Title2 = strcat('Efficiency  vs. Vertical Speed, with \eta_{hover}=',num2str(etha_emp_plot_design_Hov));
    title(Title2)
    xlabel('Vertical Speed - V_v (m/s)')
    ylabel('Efficiency \eta_{climb}')
    if MATLAB_in == 1
        legend(leg)
    else
        h_legend = legend(leg)
        set(h_legend, 'Location','Best','FontSize',LFS)
    end
    grid on

    subplot(2,2,4)
    plot(Vv_vect,(Battery_mass+Battery_mass_Hov))
    Title2 = strcat('Total Battery Mass  vs. Vertical Speed, with m_{battery,hover}=',num2str(Battery_mass_Hov),' kg for t_{hover}=',num2str(time_Hov),'seg');
    title(Title2)
    xlabel('Vertical Speed - V_v (m/s)')
    ylabel('Battery Mass (kg) - Climb and Hoover')
    if MATLAB_in == 1
        legend(leg)
    else
        h_legend = legend(leg)
        set(h_legend, 'Location','Best','FontSize',LFS)
    end
    grid on
    Fig = Fig + 1;
end