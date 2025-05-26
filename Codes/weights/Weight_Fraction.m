function Weights = Weight_Fraction(seg,Propulsion,conv_UNITS,Aero_TH,type_missions,Weight_tier,Geo_tier,AC_CONFIGURATION)

%% identifies the aerodynamic surfaces being used
case_AC = AC_CONFIGURATION.case_AC;

%% Enter the number of mission segments
% - 1 Taxy
% - 2 TakeOff
% - 3 = Climb
% - 4 = VTOL Climb
% - 5 = Cruise
% - 6 = Load Deployment
% - 7 = Turn
% - 8 = Descent
% - 9 = Descent (VTOL)
% - 10 = Alternative Airport climb to 3000ft
% - 11 = Turn loitter 45 min
% - 12 = Landing
% - 13 = Dummy to complete the 3 segment requirement for AP

CD0 = Aero_TH.CD0;
CD1 = Aero_TH.CD1;
CD2 = Aero_TH.CD2;

AR_w1 = Geo_tier.AR_w1;
S_w1 = Geo_tier.S_w1;

% Conversión desity to imperial
rho2imp = 0.06242796;
mps2ftps = conv_UNITS.mps2ftps; 
lb2kg = conv_UNITS.lb2kg;
kg2lb = conv_UNITS.kg2lb;
g = conv_UNITS.g;
mps2knot = conv_UNITS.mps2knot;
m22ft2 = conv_UNITS.m22ft2;

% Select input data
% propeller especific fuel consumptio (lb/hr/bhp)
% cbhp_1 = 0.8;
% cbhp_2 = 0.9;
% cbhp_3 = 1.0;
% cbhp_4 = 1.1;
% cbhp_5 = 1.2;
% cbhp_6 = 1.36;

switch case_AC
    case 1 % case_AC = 1 - EMERGENTIA 1:1
        % % Propulsive efficiency
        eta_mp = Propulsion.eta_mp;
        % Select propeller especific fuel consumptio (lb/hr/bhp)
        cbhp = 0.9;
        %% Empty Weight Fraction from Nicolai
        % Estimación fraciones de peso en vacio
        % propEndurance < 12 hrs
        A = 2.18;
        Cof = 0.815-1;
        % propEndurance > 12 hrs
        A = 1.66;
        Cof = 0.815-1;
        % Turbine ISR Endurance
        A = 2.78;
        Cof = 0.815-1;
        % Turbine Maneuver UCAV
        A = 3.53;
        Cof = 0.815-1;
        
        %% Empty Weight Fraction from Raymer
        % Estimación fraciones de peso en vacio
        % UAV - Tac Recce & UCAV
        A = 1.67;
        Cof = - 0.16;
        % UAV - high altitude
        A = 2.75;
        Cof = - 0.18;
        % UAV - small
        A = 0.97;
        Cof = - 0.06;
        
    case 2 % case_AC = 2 - EMERGENTIA 1:2
        % % Propulsive efficiency
        eta_mp = Propulsion.eta_mp;
        % Select propeller especific fuel consumptio (lb/hr/bhp)
        cbhp = 0.9;
        
        %% Empty Weight Fraction from Nicolai
        % Estimación fraciones de peso en vacio
        % propEndurance < 12 hrs
        A = 2.18;
        Cof = 0.815-1;
        % propEndurance > 12 hrs
        A = 1.66;
        Cof = 0.815-1;
        % Turbine ISR Endurance
        A = 2.78;
        Cof = 0.815-1;
        % Turbine Maneuver UCAV
        A = 3.53;
        Cof = 0.815-1;
        
        %% Empty Weight Fraction from Raymer
        % Estimación fraciones de peso en vacio
        % UAV - Tac Recce & UCAV
        A = 1.67;
        Cof = - 0.16;
        
        % UAV - high altitude
        A = 2.75;
        Cof = - 0.18;
        
        % UAV - small
        A = 0.97;
        Cof = - 0.06;
    case 3 % case_AC = 3 - PEPIÑO XXL
        % % Propulsive efficiency
        eta_mp = Propulsion.eta_mp;
        % Select propeller especific fuel consumptio (lb/hr/bhp)
        cbhp = 0.9;
               %% Empty Weight Fraction from Nicolai
        % Estimación fraciones de peso en vacio
        % propEndurance < 12 hrs
        A = 2.18;
        Cof = 0.815-1;
        % propEndurance > 12 hrs
        A = 1.66;
        Cof = 0.815-1;
        % Turbine ISR Endurance
        A = 2.78;
        Cof = 0.815-1;
        % Turbine Maneuver UCAV
        A = 3.53;
        Cof = 0.815-1;
        
        %% Empty Weight Fraction from Raymer
        % Estimación fraciones de peso en vacio
        % UAV - Tac Recce & UCAV
        A = 1.67;
        Cof = - 0.16;
        
        % UAV - high altitude
        A = 2.75;
        Cof = - 0.18;
        
        % UAV - small
        A = 0.97;
        Cof = - 0.06;
    case 4 % case_AC = 4 - Comercial
        % % Propulsive efficiency
        eta_mp = Propulsion.eta_mp;
        % Select propeller especific fuel consumptio (lb/hr/bhp)
        cbhp = 0.9;
        %% Empty Weight Fraction from Raymer
        % Twin Turboprop
        A = 0.97;
        Cof = - 0.05;
    case 5 % case_AC = 5 - WIG
        % % Propulsive efficiency
        eta_mp = 0.7;
        % Select propeller especific fuel consumptio (lb/hr/bhp)
        cbhp = 0.9;
        % Flying Boat
        A = 1.09;
        Cof = - 0.05;
        
        a_coef = 0;
        b_coef  = 0.42;
        C1 = -0.01;
        C2 = 0.10;
        C3 = 0.05;
        C4 = -0.12;
        C5 = 0.18;
        
        hp = 8500*2;
end

% Flight level data
N = length(type_missions);

%____Cálculo de las fracciones de peso para calcular el W0
% Cálculos de la fración wx/w0
% La misión consiste en 5 segments
for i=1:N
    if seg(i).datos.mision == 1; % Taxy
        str_mission{i} = strcat('Taxy/ ');
        wiw1m(i) = 0.97;
    elseif seg(i).datos.mision == 2; % Despegue
        str_mission{i} = strcat('Despegue/ ');
        wiw1m(i) = 0.97;
    elseif seg(i).datos.mision == 3; % Subida
        str_mission{i} = strcat('Subida/ ');
        M = seg(i).datos.Mach;
%         wiw1m(i) = 1.0065 - 0.0325*M;
        wiw1m(i) = 0.985;
        %     if M > 0.1
        %         segpes2 = 1.0065 - 0.0325*M;
        %     else
        %         segmes2 = 1;
        %     end
    elseif seg(i).datos.mision == 5; % Crucero Range
        str_mission{i} = strcat('Crucero Range/ ');
        R = seg(i).datos.dist_final*mps2ftps;
        M = seg(i).datos.Mach;
        a = seg(i).datos.Data_ATM.a;
        V = M*a;
        V_knot = V*mps2knot;
        V_maxknot = 1.25*V_knot;
        rho = seg(i).datos.Data_ATM.rho;
        Vfps = V*mps2ftps;
        C = cbhp*(1/(60*60))*Vfps/(550*eta_mp);
        %     q = 0.5*rho*rho2imp*Vfps^2;
        %     L_D = 1/((q*CD0/(W_S)) + W_S*CD2/q);
        % Optimo
        L_D = 1/(2*sqrt(CD0*CD2)); % Max Range Prop
        %     L_D = 0.943*(1/(2*sqrt(CD0*CD2))); % Max Range Jet
        wiw1m(i) = exp(-R*(1/((Vfps/C)*L_D)));
        
        % Estimation Wing Loading for Cruise
        %------------------- Cálculo de la Carga Alar
        q = 0.5*rho*(V^2);
        CL_Opt_Rng=sqrt(CD0/CD2); % max Range prop
        W_S_Cruise=(q*CL_Opt_Rng)/g; %[kg/m^2]
                        
        wxw0_im2 = 1;
        for j=1:i
            wxw0 = wiw1m(j)*wxw0_im2;
            wxw0_im2 = wxw0;
        end
        WXW0 = wxw0_im2;
        W_S_Cruise_cr = W_S_Cruise*(1/WXW0);
                
    elseif seg(i).datos.mision == 7; % Turn - Endurance
        str_mission{i} = strcat('Crucero Endurance/ ');
        n_lt = 1.00;
        E = seg(i).datos.t_final;
        M = seg(i).datos.Mach;
        a = seg(i).datos.Data_ATM.a;
        V = M*a;
        Vfps = V*mps2ftps;
        C = cbhp*(1/(60*60))*Vfps/(550*eta_mp);
        %     q = 0.5*rho*Vfps^2;
        %     L_D = 1/((q*CD0/(W_S)) + W_S*CD2/q);
        % Optimo
        L_D = 0.866*(1/(2*sqrt(CD0*CD2))); % Max Endurance Prop
        %     L_D = (1/(2*sqrt(CD0*CD2))); % Max Endurance Jet
        % Constant Velocity and constant CL
        temp = exp(-E*(1/((1/C)*L_D)));
%         w1w1m(i) = exp(-E*(1/((1/C)*L_D)))
        wiw1m(i) = temp;
        % Constant Velocity and constant altitude
        %     w1w1m = 1-(2/(1+(1/tan(E*(1/((1/(C*n_lt))*(1/sqrt(CD0*CD2))))))));
        
                % Estimation Wing Loading for Cruise
        %------------------- Cálculo de la Carga Alar
        q = 0.5*rho*(V^2);
        CL_Opt_End=sqrt(3*CD0/CD2); % max Range prop
        W_S_Endurance=(q*CL_Opt_End)/g; %[kg/m^2]
                        
        wxw0_im2 = 1;
        for j=1:i
            wxw0 = wiw1m(j)*wxw0_im2;
            wxw0_im2 = wxw0;
        end
        WXW0 = wxw0_im2;
        W_S_Endurance_cr = W_S_Endurance*(1/WXW0);

    elseif seg(i).datos.mision == 8; % Descenso
        str_mission{i} = strcat('Descenso/ ');
        wiw1m(i) = 0.993;
        
    elseif seg(i).datos.mision == 12; % Aterrizaje
        str_mission{i} = strcat('Aterrizaje/ ');
        wiw1m(i) = 0.995;
    end
    W1W1M(i) = wiw1m(i);
end

% fprintf('Selected Mission Combination')
% str_mission
% for n=1:
%     STR_mission(n) = str_mission{i};
% end
% fprintf(STR_mission)

% Relación de pesos de misión
wxw0_im1 = 1;
for j=1:N
    wxw0 = wiw1m(j)*wxw0_im1;
    wxw0_im1 = wxw0;
end
WxW0 = wxw0_im1;
    
% Relación fuel-peso inicial
WfW0 = 1.06*(1-WxW0); % reserva 6%

m_payload = Weight_tier.m_payload;
Wpayload = m_payload*kg2lb; % lb; 

% fprintf('-------------------')

Delta_W_error = 1;
W_error = 10;
W0_guess = 1000; % lbs 
count=0;
while abs(W_error) > Delta_W_error
    if count < 1
%         hp_W0_guess = 0.03*V_maxknot^0.23;
%         WeW0 = a_coef + b_coef*(W0_guess^(C1))*(AR_w1^C2)*((hp_W0_guess)^C3)*((W0_guess/(S_w1*m22ft2))^C4)*(V_maxknot^C5); % lb
        %         WeW0 = A*W0_guess^(Cof) % lb
        WeW0 = 10^(-0.0920)*W0_guess^(0.9713-1);
    else
%         WeW0 = a_coef + b_coef*(W0_guess^(C1))*(AR_w1^C2)*((hp_W0_guess)^C3)*((W0_guess/(S_w1*m22ft2))^C4)*(V_maxknot^C5); % lb
        %         WeW0 = A*W0_guess^(Cof) % lb
        WeW0 = 10^(-0.0920)*W0_guess^(0.9713-1); % lb
    end
%     WeW0
%     WfW0
    WpW0 = 1 - (WeW0) - WfW0;
    W0 = Wpayload / WpW0;
    W0_guess = (W0 + W0_guess)/2;
    W_error = W0 - W0_guess;
    count = count+1;
end

m_T0 = W0*lb2kg;
m_f = WfW0*W0*lb2kg;
m_e = WeW0*W0*lb2kg;

% mass
Weights.m_T0 = m_T0;
Weights.m_f = m_f;
Weights.m_e = m_e;
Weights.m_p = m_payload;
Weights.m_f_W0 = m_f/m_T0;
Weights.m_e_W0 = m_e/m_T0;

%------------------- Cálculo de los Cl óptimos de crucero y Loiter
 
% CL_Opt_End=sqrt(3*CD0/CD2)
% 
% CL_Opt_Rng=sqrt(CD0/CD2)
%  
% %------------------- Cálculo de la Carga Alar
    % q = 0.5*rho*(v^2);
% W_S_Cruise=(0.5*rho*(v^2)*CL_Opt_Rng)/g %[kg/m^2]
%  
% V_loi=v/1.2 % aproximación v_crucero aprox 1.2 V loiter
%  
% W_S_End=(0.5*ro*(V_loi^2)*CL_Opt_End)/g  %[kg/m^2]
 
% Tomamos la menor carga alar
 
% if W_S_Cruise< W_S_End
%     W_S=W_S_Cruise;
% else
%     W_S=W_S_End;
% end
% W_S
 
%-------------------- Características del Ala
 
S_w1=(m_T0)/W_S_Cruise_cr; % m^2
 
AR_w1 = Geo_tier.AR_w1;
b_w1=sqrt(AR_w1*S_w1); %[m]
 
c_w1 = S_w1/b_w1; %[m]    

% ----------------- Características del Perfil elegido
% 
% % Perfil ST. CYR 24 (Bartel)
% 
% % Características 2D
% 
% cl_0=0.8;
% cl_max=1.65
% 
% % Características 2D con 30º de Flap
% 
% cl_0_flap=1.55
% cl_max_flap=2.09
% 
% % Cl del ala con configuración limpia (aproximación)
% CL_ala=0.9*cl_0
% CL_Max_ala=0.9*cl_max

%---------------- Velocidades Características
% 
% fprintf('velocidades Características')
% 
% K=1/(pi*AR*e);
% 
% fprintf('H=4000m')
% H=4000;
% 
% [Temp,A,P,RHO]=atmosisa(H);
% 
% V_Stall=sqrt((2*W0final*g)/(S*RHO*CL_Max_ala))
% 
% V_Min=1.2*V_Stall
% 
% V_min_pot=sqrt((2*W0final*g)/(RHO*S)*sqrt(K/(3*Cd0)))
% 
% V_max_range=sqrt((2*W0final*g)/(RHO*S)*sqrt(K/(Cd0)))
% 
% fprintf('H=0m')
% H=0;
% 
% [Temp,A,P,RHO]=atmosisa(H);
% 
% V_Stall=sqrt((2*W0final*g)/(S*RHO*CL_Max_ala))
% 
% V_Min=1.2*V_Stall
% 
% V_min_pot=sqrt((2*W0final*g)/(RHO*S)*sqrt(K/(3*Cd0)))
% 
% V_max_range=sqrt((2*W0final*g)/(RHO*S)*sqrt(K/(Cd0)))