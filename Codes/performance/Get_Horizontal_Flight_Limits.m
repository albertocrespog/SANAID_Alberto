function [Fig] = Get_Horizontal_Flight_Limits(Aero_TH,pp,RPMMAX_APC,Prop_data,Weight_tier,conv_UNITS,Performance_preliminar,...
    Geo_tier,AC_CONFIGURATION,N_Vh,Vh_vect,Fig,PLOT_Get_Horizontal_Flight_Limits,solve_w_fero,Plot_Options)        

rho = Performance_preliminar.rho;
S_ref = Geo_tier.S_ref;
dist = Performance_preliminar.Range;
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

S_ref = Geo_tier.S_ref;
g = conv_UNITS.g;

% Initialize Figs
Fig = Fig + 1;

for i = 1:N_Vh
    Vh = Vh_vect(i);
     % Estimation of Desired Thrust
    q_inf = 0.5*rho*Vh^2;
    m_TOW = Weight_tier.m_TOW;
    
    q_inf = 0.5*rho*Vh^2;
    CL(i) = (m_TOW*g)/(q_inf*S_ref);
    CD0 = Aero_TH.CD0;
    CD1 = Aero_TH.CD1;
    CD2 = Aero_TH.CD2;
    CD = CD0 + CD1*CL(i) + CD2*CL(i)^2;
    Dh = CD*q_inf*S_ref;
    
    D_prop_h(i) = determination_D_prop_horizontal(Aero_TH,pp,RPMMAX_APC,Vh,Prop_data,Weight_tier,conv_UNITS,Performance_preliminar,Geo_tier,AC_CONFIGURATION,Dh);
    D_prop_inch_h(i) = D_prop_h(i)*100/2.54;
    

    Fdes = Dh;
    n(i) = Selection_J_CT_F_design(Vh,Prop_data,Fdes,rho,D_prop_h(i),solve_w_fero);
    RPM(i) = n(i)*60;
    
    Propulsion_h{i} = get_EngineProperties_design(Vh,rho,n(i),Prop_data,AC_CONFIGURATION,D_prop_h(i));
    
    RPM_plot_design_h(i) = Propulsion_h{i}.RPM;
    T_eng_plot_design_h(i) = Propulsion_h{i}.Ti_eng;
    Pe_eng_plot_design_h(i) = Propulsion_h{i}.Pe_eng;
    Pe_plot_design(i) = Propulsion_h{i}.Pe;
    etha_ep_plot_design(i) = Propulsion_h{i}.etha_ep;
    etha_em_plot_design(i) = Propulsion_h{i}.etha_em;
    etha_emp_plot_design(i) = Propulsion_h{i}.etha_emp;
    time_H(i) = dist/Vh;
    Energy_eng_plot_design_h(i) = (1/(1-tau_reserve))*Propulsion_h{i}.Pe*time_H(i)/3600; % En W-h
    Energy_eng_plot_design_h_kJ(i) = (1/(1-tau_reserve))*Propulsion_h{i}.Pe*time_H(i)/1000; % En kJ
    Battery_mass(i) = Energy_eng_plot_design_h(i)/(SE_battery);
end

LS = Plot_Options.LS;
FS = Plot_Options.FS;
LFS = Plot_Options.LFS*0.75;
leg = strcat('W_{payload}= ',num2str(Weight_tier.M_PAYLOAD.m_cargo),'kg & Safety Factor =',num2str(pp));

if PLOT_Get_Horizontal_Flight_Limits == 1
    
    figure(Fig)
    subplot(2,2,1)
    plot(Vh_vect,D_prop_inch_h)
    title('Prop Diameter vs. Horizontal Speed')
    xlabel('Horizontal Speed - V_h (m/s)')
    ylabel('Prop Diameter - D (inch)')
    if MATLAB_in == 1
        legend(leg)
    else
        h_legend = legend(leg)
        set(h_legend, 'Location','Best','FontSize',LFS)
    end
    grid on
    
    subplot(2,2,2)
    plot(Vh_vect,T_eng_plot_design_h)
    title('Thrust per Engine vs. Horizontal Speed')
    xlabel('Horizontal Speed - V_h (m/s)')
    ylabel('Thrust per Engine - (T_{eng}) (N)')
    if MATLAB_in == 1
        legend(leg)
    else
        h_legend = legend(leg)
        set(h_legend, 'Location','Best','FontSize',LFS)
    end
    grid on
    
    subplot(2,2,3)
    plot(Vh_vect,Pe_plot_design/1000)
    title('Total Electric Power vs. Horizontal Speed')
    xlabel('Horizontal Speed - V_v (m/s)')
    ylabel('Total Electric Power - (P_{eng})(KW)')
    if MATLAB_in == 1
        legend(leg)
    else
        h_legend = legend(leg)
        set(h_legend, 'Location','Best','FontSize',LFS)
    end
    grid on
    
    subplot(2,2,4)
    plot(Vh_vect,Pe_eng_plot_design_h/1000)
    title('Electric Power per Engine vs. Horizontal Speed')
    xlabel('Horizontal Speed - V_h (m/s)')
    ylabel('Electric Power per Engine - (P_{eng})(KW)')
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
    plot(Vh_vect,Energy_eng_plot_design_h/1000)
    title('Energy vs. Horozontal Speed')
    xlabel('Horizontal Speed - V_h (m/s)')
    ylabel('Energy - E_{h} (kW-h)')
    if MATLAB_in == 1
        legend(leg)
    else
        h_legend = legend(leg)
        set(h_legend, 'Location','Best','FontSize',LFS)
    end
    grid on

    subplot(2,2,2)
    plot(Vh_vect,etha_emp_plot_design)
    title('Efficiency  vs. Horizontal Speed')
    xlabel('Horizontal Speed - V_h (m/s)')
    ylabel('Efficiency \eta_{h}')
    if MATLAB_in == 1
        legend(leg)
    else
        h_legend = legend(leg)
        set(h_legend, 'Location','Best','FontSize',LFS)
    end
    grid on

    subplot(2,2,3)
    plot(Vh_vect,(Battery_mass))
    title('Battery Mass vs. Horizontal Speed')
    xlabel('Horizontal Speed - V_h (m/s)')
    ylabel('Battery Mass (kg)')
    if MATLAB_in == 1
        legend(leg)
    else
        h_legend = legend(leg)
        set(h_legend, 'Location','Best','FontSize',LFS)
    end
    grid on

    subplot(2,2,4)
    plot(Vh_vect,RPM_plot_design_h)
    title('RPMs vs Horizontal Speed')
    xlabel('Horizontal Speed - V_h (m/s)')
    ylabel('RPMs')
    grid on
    
    Fig = Fig + 1;
    figure(Fig)
    plot(Vh_vect,CL)
    title('Lift Coefficient vs. Horizontal Speed')
    xlabel('Horizontal Speed - V_h (m/s)')
    ylabel('C_L')
    grid on
    Fig = Fig + 1;
end
