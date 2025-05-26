function [T,c,Vnmax,delta_act] = propulsive_model_v5(V,h,delta_p)

%Mediante esta rutina se calcula la tracción y el consumo para una
%determinada altitud de vuelo, posición de palanca de gases y diámetro de
%la hélice.
%
%La planta motora considerada es un motor Hirth 4103, con un modelo de
%potencia cúbico, donde el único parámetro estimado es nPzero.
%
%Se considera un modelo de consumo específico lineal en las revoluciones.
%
%Última actualización: 18/02/2016 (Se ha cambiado la función de splines).

%Determinación de la penalización por altitud: Modelo de Raymer.
delta_rho = (1-6.5e-3*h/288.15)^4.25193; %Troposfera ISA.
rho       = 1.225*delta_rho;
delta_h   = delta_rho-(1-delta_rho)/7.55;
%Resultados en consonancia con lo que establece Rotax para su modelo 912
%iS de inyección electrónica (ver manual de operación

%Datos de la planta propulsora: Hélice 5868-9, Clark Y sections, 2 blades,
%beta_n = 15º.
Jdata    = [(0:0.1:0.7)'; 0.75; 0.8; 0.84];
cPdata   = [0.041; 0.041; 0.041; 0.04; 0.037; 0.032; 0.025; 0.016; 0.0106; 0.0048; 0];
cPspline = spline(Jdata(3:end)',cPdata(3:end)')';
cTdata   = [0.099; 0.096; 0.089; 0.078; 0.064; 0.049; 0.034; 0.0185; 0.0106; 0.0025; -0.0042];
cTspline = spline(Jdata',cTdata')';
%Diámetro de hélice escogido
D        = 28*0.0254; %[m]

%Corrección de John T. Lowry (Slow Down Efficiency Factor):
Z    = 2*sqrt(0.122/pi)/D; %Diámetro fuselaje / diámetro hélice.
SDEF = 1.05263-0.04185*Z-0.01481*Z^2-0.62001*Z^3;

%Datos de la planta motora: Hirth 4103/4102 - Nuevo manual.
Pmot   = 5740; 
nPzero = 30;           nPmax  = 116;
nmin   = 30;           nmax   = 6500/60;

%En realidad es el starter/generator el que marca las revoluciones máximas
%de modo que nmax no es 6500/60, sino:
nmax   = 6059/60;

%Pmax   = delta_p*delta_h*Pmot;
Pcoef  = polyfit([nPzero; 79; nPmax],delta_p*delta_h*[0; 3730; Pmot],2);

%Datos del alternador/arrancador del Hirth:
dPedn  = 500/50.49;     %A 50.49 rev/s produce 500 W.
Pelec  = 1.1*115;       %Estimación mayorada de consumos eléctricos.
Pcomax = Pelec/0.75;    %Se incluye rendimiento de transformaciones.
%nPcons = Pelec/dPedn;   %Mínima velocidad de giro a la que se saca Pelec.
nPcons = 20;   %Mínima velocidad de giro a la que se saca Pelec.

%Determinación del quiebro de la curva (punto break)
Pnmax    = Pcoef*[nmax^2; nmax; 1];
cPnmax   = (Pnmax-Pcomax)/(rho*nmax^3*D^5);
Jnmax    = spline(cPdata(4:end),Jdata(4:end),cPnmax);
Vnmax    = Jnmax*nmax*D;
misops   = optimset('Display','off','TolX',1e-9);

%Determinación de cómo está operando el motor
if V<Vnmax,
    %Operación a delta_act = delta_p:
    delta_act = delta_p;

    Jguess = V/(nmax*D);
    [J,fval,flag1] = fzero(@(x) equilfun(x,V,Jdata,cPdata,cPspline,Pcoef,nPcons,Pcomax,rho,D),Jguess,misops);
    if (flag1>0)&&(V/(J*D)>=nPcons),
        %Se extrae la potencia necesaria para que los sistemas eléctricos
        %reciban su potencia nominal, y el método converge a una solución
        %válida:
        n    = V/(J*D);
        cT   = fnval(cTspline,J);
        T    = SDEF*cT*rho*n^2*D^4;
        Ptot = Pcoef*[n^2; n; 1];
    else
        %No se puede extraer la potencia necesaria para que los sistemas
        %eléctricos reciban su potencia nominal, sino que se extrae menos.
        Jguess = V/(nPcons*D);
        [J,fval,flag2] = fzero(@(x) equilfun(x,V,Jdata,cPdata,cPspline,Pcoef,nPcons,Pcomax,rho,D),Jguess,misops);
        if (flag2>0)&&(V/(J*D)>=nmin),
            %El método converge a una solución válida:
            n    = V/(J*D);
            cT   = fnval(cTspline,J);
            T    = SDEF*cT*rho*n^2*D^4;
            Ptot = Pcoef*[n^2; n; 1];
        else
            %En realidad el motor no da potencia, y el único punto de equilibrio posible es J2.
            J    = Jdata(end);
            n    = V/(J*D);
            cT   = cTdata(end);
            T    = SDEF*cT*rho*n^2*D^4;
            Ptot = 0;
        end
    end
else
    %Operación a n=nmax;
    n         = nmax;
    J         = V/(n*D);
    cT        = fnval(cTspline,J);
    T         = SDEF*cT*rho*n^2*D^4;
    cP        = fnval(cPspline,J);
    Ptot      = rho*nmax^3*D^5*cP + Pcomax;
    P_ratio   = Ptot/Pnmax;
    delta_act = P_ratio*delta_p;
end 

%DETERMINACIÓN DEL CONSUMO DE COMBUSTIBLE:
% La ley de consumo específico se da como ajuste lineal a partir de los
% datos del Hirth F-36, de 3000 a 6000 rpm:
% nmot_vec = (3000:250:6000)'/60;
% Pmot_vec = [6; 6.8; 8.6; 9.4; 9.7; 11; 13; 13.4; 13.7; 13.8; 14; 14.2; 15]*(6000/15);
% c_dat = [0.89; 0.99; 1.08; 1.2; 1.48; 1.6; 1.77; 1.87; 2; 2; 2.27; 2.27; 2.4];
% c_sca = (6000/(15*745.699872))*3.78541178*0.75/3600*c_dat;
% %         escala*(litros/gallon)*(kg/litro)*(hora/segundos)*c_dat;
% cE_vec = c_sca./Pmot_vec;
% cE_coef = [ones(size(cE_vec)),nmot_vec]\cE_vec;
% 
% c    = 1.2*(250/285)*(cE_coef(1)+cE_coef(2)*n)*Ptot;
%Se añade un factor de mejora tecnológica al pasar de un motor con
%carburador a un motor con inyección electrónica (Fuente: Rotax 912 S y
%Rotax 912 iS), y se añade un factor de seguridad de 1.2 para absorber la
%incertidumbre en los datos del motor y de la hélice.

%Valor proporcionado por AERTEC:
c    = 1.1*0.60/1e3/3600*Ptot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function elerr = equilfun(J,V,Jdata,cPdata,cPspline,Pcoef,nPcons,Pcomax,rho,D)

if J<=Jdata(3),
    cP = cPdata(1);
else
    cP = fnval(cPspline,J);
end

if V/(J*D)>nPcons,
    P = Pcoef*[(V/(J*D))^2; V/(J*D); 1]-Pcomax;
else
    P = Pcoef*[(V/(J*D))^2; V/(J*D); 1]-Pcomax*V/(J*D*nPcons);
end

elerr  = rho*D^2*V^3*cP - P*J^3;
