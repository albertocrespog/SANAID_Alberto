function [Propulsion] =get_EngineProperties_AP(h,V,alpha,beta,Geo_tier,Fdes,AC_CONFIGURATION)

%% Propulsion information
% 1: TIPO DE MOTOR --> 1_TURBOFAN 2_TURBOPROP 3_PISTON 4_ELECTRIC_PROP
% 2: NUMERO DE MOTORES
% 3: EMPUJE/POTENCIA A NIVEL DEL MAR: % Thrust (lbf) or Power (shp) per engine
% 4: CONSUMO ESPECIFICO: % Specificl Fuel consumption  (lb/lbf*h)/(lb/shp*h)
% 5: DERIVACION(TURBOFANES) By-pass
% 6: EFICIENCIA DE LA HELICE (ETA_P)
% 7: AVION CIVIL = 1/MILITAR = 2 - Normativa
% 8: CAPACIDAD CALORIFICA DEL COMBUSTIBLE
% 9: DIAMETRO DEL MOTOR(TURBOHELICES)

h_SL = 0;
[Temp_SL,rho_SL,p_SL,a_SL]=atmos_inter_mio(h_SL);

[Temp,rho,p,a]=atmos_inter_mio(h);

S_ref = Geo_tier.S_ref;
cmac = Geo_tier.cmac_w1;
qbar = 0.5*rho*V^2;

propul = AC_CONFIGURATION.propulsion;

% propul(1) = 3; % Type of engine
% propul(2) = n_eng; % Number of engines
% propul(3) = 5; % Thrust (lbf) or Power (shp) per engine
% propul(4) = 0.901; % Specificl Fuel consumption  (lb/lbf*h)/(lb/shp*h)
% propul(5) = 0; % By-pass
% propul(6) = 0.82; % Prop efficiency
% propul(7) = 1; % Normativa
% propul(8) = 0; % Capacidad calorifica del combustible
% propul(9) = 28*2.54/100; % Diámetro de la hélice

% PROPULSION
tipo_motor = propul(1);
n_eng = propul(2);
T_SL = propul(3);
c_SL = propul(4);
civil_mil = propul(5);
eta_p = propul(6);
derivacion = propul(7);

% Coefficients correctores consumo específico
a_1=3.559957437510763;
a_2=-10.739698199171459;
a_3= 11.989635150373475;
a_4=-5.869876557884609;
a_5=2.059994459180667;

% Calculates the delta as a function of desired Thrust (Drag)
switch tipo_motor
    case 1
        delta_T = (Fdes/n_eng)/(T_SL*(1+0.2*M^2)^(1.4./0.4)*(1-0.49.*sqrt(M))*(rho/rho_SL));
        correccion = (a_1.*delta_T.^4+a_2.*delta_T.^3+a_3.*delta_T.^2+a_4.*delta_T+a_5);
        switch derivacion
            case 1
                C =  @(W) c_SL.*correccion.*(1+1.2*M).*sqrt(Temp/Temp_SL);
            case 2
                C =  @(W) c_SL.*correccion.*(1+0.33*M).*sqrt(Temp/Temp_SL);
            case 3
                C =  @(W) c_SL.*correccion.*(1+0.16875*M).*sqrt(Temp/Temp_SL);
        end
    case 2
        delta_T = (Fdes/n_eng)/(T_SL*eta_p/V*(1+0.2.*M^2)^(1.4./0.4)*(rho*Temp/(rho_SL*Temp_SL)));
        correccion = (a_1.*delta_T.^4+a_2.*delta_T.^3+a_3.*delta_T.^2+a_4.*delta_T+a_5);
        C = c_SL.*correccion.*(1+1.44*M).*sqrt(Temp/Temp_SL).* V./eta_p;
    case 3
        delta_T = (Fdes/n_eng)/(T_SL*eta_p/V*((8.55*rho/rho_SL - 1)/7.55));
        correccion = (a_1.*delta_T.^4+a_2.*delta_T.^3+a_3.*delta_T.^2+a_4.*delta_T+a_5);
        C = c_SL.*correccion.*V./eta_p;
end
Propulsion.delta_T = delta_T;
Propulsion.C = C;
Propulsion.correccion = correccion;

switch tipo_motor
    case 1
        Ti_eng = delta_T.*(T_SL.*(1+0.2.*M.^2).^(1.4./0.4).*(1-0.49.*sqrt(M)).*(rho./rho_SL));
        Ti = Ti_eng*n_eng;
        CT = T / (qbar*Sref);
        dT_dV = 1.400000000*T_SL*(1+.2*V^2/a^2)^2.500000000*(1-.49*sqrt(V/a))*rho*V/(rho_SL*a^2)-...
            .2450000000*T_SL*(1+.2*V^2/a^2)^3.500000000*rho/(sqrt(V/a)*a*rho_SL);
        dCT_dV = 2.800000000*T_SL*(1+.2*V^2/a^2)^2.500000000*(1-.49*sqrt(V/a))/(rho_SL*V*S_ref*a^2)-...
            .4900000000*T_SL*(1+.2*V^2/a^2)^3.500000000/(sqrt(V/a)*a*rho_SL*V^2*S_ref)-...
            4.000000000*T_SL*(1+.2*V^2/a^2)^3.500000000*(1-.49*sqrt(V/a))/(rho_SL*V^3*S_ref);
    case 2
        Ti_eng = delta_T.*(T_SL.*eta_p./V.*(1+0.2.*M.^2).^(1.4./0.4).*(rho.*Temp./(rho_SL*Temp_SL)));
        Ti = Ti_eng*n_eng;
        Pi_eng = Ti_eng*V;
        Pi = Pi_eng*n_eng;
        Pe = Pi/eta_p;
        Pe_eng = Pe/n_eng;
        CT = Ti_eng / (qbar*S_ref);
        CP = Pi_eng / (qbar*S_ref*cmac);
        CQ = CP/(2*pi);
        Q_eng = CQ 
        dT_dV = 1.400000000*delta_T*T_SL*eta_p*(1+.2*V^2/a^2)^2.500000000*rho*Temp/(rho_SL*Temp_SL*a^2)-...
            delta_T*T_SL*eta_p*(1+.2*V^2/a^2)^3.500000000*rho*Temp/(V^2*rho_SL*Temp_SL);
        dCT_dV = 2.800000000*delta_T*T_SL*eta_p*(1+.2*V^2/a^2)^2.500000000*Temp/(V^2*rho_SL*Temp_SL*S_ref*a^2)-...
            6.000000000*delta_T*T_SL*eta_p*(1+.2*V^2/a^2)^3.500000000*Temp/(V^4*rho_SL*Temp_SL*S_ref);
        dP_dV = 1.400000000*delta_T*T_SL*eta_p*(1+.2*V^2/a^2)^2.500000000*rho*Temp*V/(rho_SL*Temp_SL*a^2);
        dCP_dV = 2.800000000*delta_T*T_SL*eta_p*(1+.2*V^2/a^2)^2.500000000*Temp/(rho_SL*Temp_SL*V*S_ref*cmac*a^2)-...
            4.000000000*delta_T*T_SL*eta_p*(1+.2*V^2/a^2)^3.500000000*Temp/(rho_SL*Temp_SL*V^3*S_ref*cmac);
        %Storing data
        Propulsion.Pi_eng = Pi_eng;
        Propulsion.Pi = Pi;
        Propulsion.Pe = Pe;
        Propulsion.Pe_eng = Pe_eng;
        Propulsion.d_P_d_V = dP_dV;
        Propulsion.d_CP_d_V = dCP_dV;
    case 3
        Ti_eng = delta_T.*(T_SL.*eta_p./V.*((8.55.*rho./rho_SL - 1)./7.55));
        Ti = Ti_eng*n_eng;
        Pi_eng = Ti_eng*V;
        Pi = Pi_eng*n_eng;
        Pe = Pi/eta_p;
        Pe_eng = Pe/n_eng;
        CT = Ti_eng / (qbar*S_ref);
        CP = Pi_eng / (qbar*S_ref*cmac);
        dT_dV = -.1324503311*delta_T*T_SL*eta_p*(8.55*rho/rho_SL-1)/V^2;
        dCT_dV = -.7947019866*delta_T*T_SL*eta_p*(8.55*rho/rho_SL-1)/(V^4*rho*S_ref);
        dP_dV = 0;
        dCP_dV = -.5298013244*delta_T*T_SL*eta_p*(8.55*rho/rho_SL-1)/(rho*V^3*S_ref*cmac);
        %Storing data
        Propulsion.Pi_eng = Pi_eng;
        Propulsion.Pi = Pi;
        Propulsion.Pe = Pe;
        Propulsion.Pe_eng = Pe_eng;
        Propulsion.d_P_d_V = dP_dV;
        Propulsion.d_CP_d_V = dCP_dV;
end
% Storing Data
Propulsion.Ti = Ti;
Propulsion.Ti_eng = Ti_eng;
Propulsion.d_T_d_V = dT_dV;
Propulsion.d_CT_d_V = dCT_dV;
