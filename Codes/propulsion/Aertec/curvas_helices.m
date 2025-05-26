% Curvas C_P y C_T de las distintas hélices

%% CURVAS CP
clear all;  clc

% HÉLICE BIPALA 36x12
CP0_3612=0.0189 ;
CP1_3612=0.01604 ;
CP2_3612=-0.06985 ;

fCP_3612=@(x) CP0_3612 + CP1_3612*x + CP2_3612*x^2;

% HÉLICE TRIPALA 28x12
CP0_2812=0.03036;
CP1_2812=0.03298 ;
CP2_2812=-0.1007 ;

fCP_2812=@(x) CP0_2812 + CP1_2812*x + CP2_2812*x^2;

% HÉLICE BIPALA 34x12
CP0_3412=0.0161;
CP1_3412=0.02106;
CP2_3412=-0.06772;

fCP_3412=@(x) CP0_3412 + CP1_3412*x + CP2_3412*x^2;

% HÉLICE TRIPALA 31x12
CP0_3112=0.02014;
CP1_3112=0.02014;
CP2_3112=-0.07632;

fCP_3112=@(x) CP0_3112 + CP1_3112*x + CP2_3112*x^2;


% FIGURA DE COMPARACIÓN DE TODAS
figure (1)

fplot(fCP_3612,[0,0.8],'b')
hold on
fplot(fCP_2812,[0,0.8],'g')
hold on
fplot(fCP_3412,[0,0.8],'r')
hold on
fplot(fCP_3112,[0,0.8],'k')
hold on
axis([0 0.8 0 0.035])
xlabel('J')
ylabel('C_P')

legend('36x12','28x12','34x12','31x12')



%% CURVAS CT
clear all;  clc

% HÉLICE BIPALA 36x12
CT0_3612=0.07897 ;
CT1_3612=-0.1111  ; 
CT2_3612=-0.03181 ;

fCT_3612=@(x) CT0_3612 + CT1_3612*x + CT2_3612*x^2;

% HÉLICE TRIPALA 28x12
CT0_2812=0.1062 ;
CT1_2812=-0.1024 ;
CT2_2812=-0.07823 ;

fCT_2812=@(x) CT0_2812 + CT1_2812*x + CT2_2812*x^2;

% HÉLICE BIPALA 34x12
CT0_3412=0.07102;
CT1_3412=-0.09107;
CT2_3412=-0.03765; 

fCT_3412=@(x) CT0_3412 + CT1_3412*x + CT2_3412*x^2;

% HÉLICE TRIPALA 31x12
CT0_3112=0.0807;
CT1_3112=-0.08869;
CT2_3112=-0.05565;

fCT_3112=@(x) CT0_3112 + CT1_3112*x + CT2_3112*x^2;


% FIGURA DE COMPARACIÓN DE TODAS
figure (2)

fplot(fCT_3612,[0,0.8],'b')
hold on
fplot(fCT_2812,[0,0.8],'g')
hold on
fplot(fCT_3412,[0,0.8],'r')
hold on
fplot(fCT_3112,[0,0.8],'k')
hold on
axis([0 0.8 0 0.11])
xlabel('J')
ylabel('C_T')

legend('36x12','28x12','34x12','31x12')




























