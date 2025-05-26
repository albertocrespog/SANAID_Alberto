function [fuel_taxi] = analisis_taxi(taxi,propul,i)

% --> %% ARRANQUE DE LOS MOTORES + WARM-UP = TAXI
%       % 1: TEMPERATURA LOCAL
%       % 2: ALTURA LOCAL
%       % 3: PRESION LOCAL
%       % 4: PALANCA DE RALENTI EN TAXI = 0.05
%       % 5: VELOCIDAD A LA QUE HACE EL TAXI
%       % 6: TIEMPO DE ESPERA EN TAXI

%% PROPUL
% 1: TIPO DE MOTOR --> 1_TURBOFAN 2_TURBOPROP 3_PISTON
% 2: NUMERO DE MOTORES
% 3: EMPUJE/POTENCIA A NIVEL DEL MAR
% 4: CONSUMO ESPECIFICO
% 5: AVION CIVIL =1/MILITAR = 2
% 6: EFICIENCIA DE LA HELICE (ETA_P)
% 7: DERIVACION(TURBOFANES)

%%

R = 287.058;
a_taxi = sqrt(1.4*R*taxi(1));
tiempo_taxi = taxi(6);
                    
switch propul(1)
        case 1
            h = taxi(2);
            delta_T = taxi(4);
            delta_T0 = 1;
            M = taxi(5)/a_taxi;
            M_0 = 0.15;
            T_SL = propul(3) * propul(2);
            c_SL = propul(4);
            derivacion = propul(7);
            
           [T_T0_taxi,C_taxi,T_taxi] = modprop_fan(h,delta_T,delta_T0,M,M_0,T_SL,c_SL,derivacion);
           fuel_taxi = C_taxi*T_taxi*tiempo_taxi;
        
        case 2
            h = taxi(2);
            delta_T = taxi(4);
            delta_T0 = 1;
            M = taxi(5)/a_taxi;
            M_0 = 0.15;
            P_SL = propul(3) * propul(2);
            cp_P_SL = propul(4);
           
           [T_T0_taxi,C_taxi,P_taxi,T_taxi]=modprop_prop(h,delta_T,delta_T0,M,M_0,P_SL,cp_P_SL);
            fuel_taxi = C_taxi*T_taxi*tiempo_taxi;
         
     case 3
            h = taxi(2);
            delta_P = taxi(4);
            delta_P0 = 1;
            M = taxi(5)/a_taxi;
            M_0 = 0.15;
            B_SL = propul(3) * propul(2);
            cp_P_SL = propul(4);
           
           [T_T0_taxi,C_taxi,B_taxi,T_taxi]=modprop_piston(h,delta_P,delta_P0,M,M_0,B_SL,cp_P_SL);
            fuel_taxi = C_taxi*T_taxi*tiempo_taxi;
end

end        