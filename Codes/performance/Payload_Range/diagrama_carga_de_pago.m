%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function diagrama_carga_de_pago(PL_max,PL_min,ZFM,m_TO,mF_max,C_L_max_CR,S,C_D0_CR,C_D1_CR,C_D2_CR,h_CR,h_f)

figure, grid on, hold on
xlabel('Autonomía [h]')
ylabel('Carga de pago [kg]')
axis([0 13 PL_min PL_max+2])
h = text(10,10,['Masa al despegue ' num2str(m_TO(end)) ' kg']);
set(h,'rotation',-40);

misops_fzero  = optimset('TolX',1e-7); %Opciones para fzero
misops_fsolve = optimset('Display','off','TolX',1e-7,'TolFun',1e-7,'MaxIter',10,'DiffMinChange',1e-8,'FinDiffType','central'); %Opciones para fsolve

%El punto A se corresponde al caso de no fuel
PL_A  = PL_max;
End_A = 0;

%El punto B se corresponde al caso de alcanzar la masa máxima al despegue
%(debe coincidir con las 8 horas de vuelo)
PL_B  = PL_max;
End_B = fzero(@(End) ZFM + PL_B - combustible_analisis(m_TO(end),C_L_max_CR,S,C_D0_CR,C_D1_CR,C_D2_CR,End*3600,h_CR,h_f),[1 15],misops_fzero);

plot([End_A End_B],[PL_A PL_B],'linewidth',1)

%El punto C se corresponde con el corte entre la curva de masa máxima al
%despegue y máxima capacidad del tanque (el combustible consumido debe
%coincidir con la capacidad máxima del tanque)
flag = -1;
fval = 1;
x_0 = [PL_min 0.9*mF_max];
while flag <= 0 || norm(fval) > 1e-3 %Si la bandera es mnor que 0, entonces es que fsolve ha tenido problemas. Reanudamos la resuloción partiendo del último punto propuesto
    %También continuamos repitiendo el proceso si la norma de la
    %función objetivo es superior a 1e-3 (aproximadamente 10 gramos)
    [x_C,fval,flag] = fsolve(@(x) [(ZFM + x(1) - combustible_analisis(m_TO(end),C_L_max_CR,S,C_D0_CR,C_D1_CR,C_D2_CR,x(2)*3600,h_CR,h_f))/ZFM
        (ZFM + x(1) - combustible_analisis(ZFM + x(1) + mF_max,C_L_max_CR,S,C_D0_CR,C_D1_CR,C_D2_CR,x(2)*3600,h_CR,h_f))/ZFM],...
        x_0,misops_fsolve); %La primera estimación está basada en que se consume aproximadamente 0.9 kg por hora
    x_0 = x_C;
end
PL_C  = x_C(1);
End_C = x_C(2);

if PL_C > PL_min
    PL = linspace(PL_B,PL_C,5);
    End(1) = End_B; End(5) = End_C;
    for j = 2:length(PL)-1
        End(j) = fzero(@(End) ZFM + PL(j) - combustible_analisis(m_TO(end),C_L_max_CR,S,C_D0_CR,C_D1_CR,C_D2_CR,End*3600,h_CR,h_f),[End_B End_C],misops_fzero);
    end
else
    PL = linspace(PL_B,PL_min,5);
    End(1) = End_B;
    for j = 2:length(PL)
        End(j) = fzero(@(End) ZFM + PL(j) - combustible_analisis(m_TO(end),C_L_max_CR,S,C_D0_CR,C_D1_CR,C_D2_CR,End*3600,h_CR,h_f),[End_B End_C],misops_fzero);
    end
end
plot(End,PL,'linewidth',1)

for i = (length(m_TO)-1):-1:1
    m_TO(i)
    %Primero determinamos cuál sería la carga de pago que se correspondería
    %con autonomía 0
    PL_End0 = combustible_analisis(m_TO(i),C_L_max_CR,S,C_D0_CR,C_D1_CR,C_D2_CR,0,h_CR,h_f) - ZFM;
    
    %Si la carga de pago que se corresponde con autonomía 0 es mayor que
    %PL_max, entonces el punto B tiene el mismo significado que antes
    if PL_End0 >= PL_max
        PL_B  = PL_max;
        End_B = fzero(@(End) ZFM + PL_B - combustible_analisis(m_TO(i),C_L_max_CR,S,C_D0_CR,C_D1_CR,C_D2_CR,End*3600,h_CR,h_f),[0 8],misops_fzero);
    else
        %Si no, ahora el punto B es el punto en el que se tiene autonomía 0
        PL_B = PL_End0;
        End_B = 0;        
    end
    
    %El punto C se corresponde con el corte entre la curva de masa máxima al
    %despegue y máxima capacidad del tanque
    flag = -1;
    fval = 1;
    x_0 = x_C; %Usamos como punto de arranque el punto C obtenido en la masa anterior
    while flag <= 0 || norm(fval) > 1e-3 %Si la bandera es menor que 0, entonces es que fsolve ha tenido problemas. Reanudamos la resuloción partiendo del último punto propuesto
        %También continuamos repitiendo el proceso si la norma de la
        %función objetivo es superior a 1e-3 (aproximadamente 10 gramos)
        [x_C,fval,flag] = fsolve(@(x) [(ZFM + x(1) - combustible_analisis(m_TO(i),C_L_max_CR,S,C_D0_CR,C_D1_CR,C_D2_CR,x(2)*3600,h_CR,h_f))/ZFM
            (ZFM + x(1) - combustible_analisis(ZFM + x(1) + mF_max,C_L_max_CR,S,C_D0_CR,C_D1_CR,C_D2_CR,x(2)*3600,h_CR,h_f))/ZFM],...
            x_0,misops_fsolve); 
        x_0 = x_C;
    end
    PL_C  = x_C(1);
    End_C = x_C(2);
    
    %Si PL_C es negativo, entonces hay que cortar en el punto en el que PL = PL_min
    if PL_C < 0
        PL_C = PL_min;
        End_C = fzero(@(End) ZFM + PL_C - combustible_analisis(m_TO(i),C_L_max_CR,S,C_D0_CR,C_D1_CR,C_D2_CR,End*3600,h_CR,h_f),[End_B End_C],misops_fzero);
    end
    
    %Representamos enre PL_B y PL_C
    PL = linspace(PL_B,PL_C,5);
    End(1) = End_B; End(5) = End_C;
    for j = 2:length(PL)-1
        End(j) = fzero(@(End) ZFM + PL(j) - combustible_analisis(m_TO(i),C_L_max_CR,S,C_D0_CR,C_D1_CR,C_D2_CR,End*3600,h_CR,h_f),[End_B End_C],misops_fzero);
    end
    plot(End,PL,'linewidth',1)
end