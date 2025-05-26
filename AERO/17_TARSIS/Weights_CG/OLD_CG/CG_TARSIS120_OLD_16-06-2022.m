function OUTPUT = CG_TARSIS120(FUEL,W_PL,RACK,n_MSL,DEPOSITO,MOTOR,CATIAT120)
% Elección de la Configuración
seleccion_motor = MOTOR; % 1 para SP210, 2 para DA215
seleccion_rack =RACK; % 1 == rack, 0 == sin rack
seleccion_misiles= n_MSL; %0,1,2,3,4
seleccion_deposito= DEPOSITO; % 1 == deposito original, 2 == depósito húmedo
llenado_deposito= FUEL;  % en tanto por 1
peso_payload = W_PL;

% CATIAT120=xlsread('Tabla extraida de CATIA.xlsx','D4:M25');

% [PESO, X, Y, Z, Ix, Iy, Iz, Ixy, Ixz, Iyz]
FUS = CATIAT120(1,:);
TAIL = CATIAT120(2,:);
NSE_LGR = CATIAT120(3,:);
PPL_LGR = CATIAT120(4,:);
PPL_LGR(5:10) = PPL_LGR(5:10)*(PPL_LGR(1)/6.752);
SP210 = CATIAT120(5,:);
DA215 = CATIAT120(6,:);
WNG_LFT = CATIAT120(7,:);
WNG_RGT = WNG_LFT;
WNG_RGT(3) = -WNG_RGT(3);
WNG_RGT(10) = -WNG_RGT(10);
WNG_RGT(8)=-WNG_RGT(8);
RCK_1=CATIAT120(8,:);
RCK_0=CATIAT120(9,:);
MSL_0=CATIAT120(10,:);
MSL_1=CATIAT120(11,:);
MSL_2=CATIAT120(12,:);
MSL_3=CATIAT120(13,:);
MSL_4=CATIAT120(14,:);
PLD=CATIAT120(15,:);
PLD(1)=peso_payload;
PLD(5:10)=PLD(5:10)*(peso_payload/3.75);
SYS=CATIAT120(20,:);

W_MAZ_ELE=CATIAT120(21,1) ;
W_PEG_TOR=CATIAT120(22,1);




if seleccion_motor==1
    MPS=SP210; % SP210
elseif seleccion_motor==2
    MPS=DA215; % DA215
else
    disp('error en seleccion del motor')
end




if seleccion_rack==1
    RACK_LFT=RCK_1; % con rack
      RACK_RGT=RACK_LFT;
      RACK_RGT(3)=-RACK_RGT(3);
      RACK_RGT(10)=-RACK_RGT(10);
      RACK_RGT(8)=-RACK_RGT(8);
elseif seleccion_rack==0
    RACK_LFT=RCK_0; % sin rack
    RACK_RGT=RCK_0;
else
    disp('error en seleccion del rack')
end



if seleccion_rack==1
  if seleccion_misiles==0
      MSL=MSL_0; % sin misiles
  elseif seleccion_misiles==1
      MSL=MSL_1; % 1 misil
  elseif seleccion_misiles==2
      MSL=MSL_2; % 2 misiles
  elseif seleccion_misiles==3
      MSL=MSL_3; % 3 misiles
  elseif seleccion_misiles==4
      MSL=MSL_4; % 4 misiles
  else
      disp('error en seleccion de los misiles')
  end
else
   disp('No se pueden cargar misiles sin Rack')
end




if seleccion_deposito==1 %original/extraíble
    GAS_FULL=CATIAT120(16,:);
    GAS_EMPTY=CATIAT120(17,:);
elseif seleccion_deposito==2 %húmedo
    GAS_FULL=CATIAT120(18,:);
    GAS_EMPTY=CATIAT120(19,:);
else
    disp('error en la selección del deposito')
end


GAS_CAP=GAS_FULL-GAS_EMPTY;
GAS=GAS_EMPTY+GAS_CAP*llenado_deposito;


%% OBTENCIÓN CENTRO DE GRAVEDAD

m=1;
x=2;
y=3;
z=4;

%debug
% FUS. =FUS*0;
% GAS=GAS*0;
% PLD=PLD*0;
% RACK_RGT=RACK_RGT*0;
% WNG_RGT=WNG_RGT*0;
% NSE_LGR=WNG_RGT*0;
% PPL_LGR=WNG_RGT*0;
% SYS=WNG_RGT*0;
% MSL=WNG_RGT*0;
% W_MAZ_ELE=0;
% W_PEG_TOR=0;

sum_masas=[FUS(m)+WNG_LFT(m)+WNG_RGT(m)+TAIL(m)+NSE_LGR(m)+PPL_LGR(m)+MPS(m)+RACK_LFT(m)+RACK_RGT(m)+MSL(m)+PLD(m)+GAS(m)+SYS(m)+W_MAZ_ELE+W_PEG_TOR];

x_cg=(FUS(m)*FUS(x) + WNG_LFT(m)*WNG_LFT(x) + WNG_RGT(m)*WNG_RGT(x) + TAIL(m)*TAIL(x) + NSE_LGR(m)*NSE_LGR(x)...
    + PPL_LGR(m)*PPL_LGR(x) + MPS(m)*MPS(x) + RACK_LFT(m)*RACK_LFT(x) + RACK_RGT(m)*RACK_RGT(x) + MSL(m)*MSL(x) + PLD(m)*PLD(x) + GAS(m)*GAS(x)...
    + SYS(m)*SYS(x))/sum_masas;
y_cg=(FUS(m)*FUS(y) + WNG_LFT(m)*WNG_LFT(y) + WNG_RGT(m)*WNG_RGT(y) + TAIL(m)*TAIL(y) + NSE_LGR(m)*NSE_LGR(y)...
    + PPL_LGR(m)*PPL_LGR(y) + MPS(m)*MPS(y) + RACK_LFT(m)*RACK_LFT(y) + RACK_RGT(m)*RACK_RGT(y) + MSL(m)*MSL(y) + PLD(m)*PLD(y) + GAS(m)*GAS(y)...
    + SYS(m)*SYS(y))/sum_masas;
z_cg=(FUS(m)*FUS(z) + WNG_LFT(m)*WNG_LFT(z) + WNG_RGT(m)*WNG_RGT(z) + TAIL(m)*TAIL(z) + NSE_LGR(m)*NSE_LGR(z)...
    + PPL_LGR(m)*PPL_LGR(z) + MPS(m)*MPS(z) + RACK_LFT(m)*RACK_LFT(z) + RACK_RGT(m)*RACK_RGT(z) + MSL(m)*MSL(z) + PLD(m)*PLD(z) + GAS(m)*GAS(z)...
    + SYS(m)*SYS(z))/sum_masas;

x_morro=490.224; % distancia del mamparo al morro
CG=[x_cg + x_morro,y_cg,z_cg]/1000;

x_cg=x_cg*0.001; 
y_cg=y_cg*0.001;
z_cg=z_cg*0.001;

 %SUMA TOTAL DE INERCIA EN MORRO

I_TOTAL = FUS(5:10) + NSE_LGR(5:10) + PPL_LGR(5:10) + MPS(5:10)  + TAIL(5:10)...
   + WNG_RGT(5:10) + WNG_LFT(5:10) + RACK_LFT(5:10) + RACK_RGT(5:10) + MSL(5:10) + PLD(5:10) + GAS(5:10) + SYS(5:10);
% STEINER - INERCIA EN CG
I_TOTALCG(1) = I_TOTAL(1) - sum_masas*(y_cg^2+z_cg^2);
I_TOTALCG(2)= I_TOTAL(2)  - sum_masas*(x_cg^2+z_cg^2);
I_TOTALCG(3) = I_TOTAL(3)  - sum_masas*(x_cg^2+y_cg^2);
I_TOTALCG(4) = I_TOTAL(4)  + sum_masas*x_cg*y_cg;
I_TOTALCG(5) = I_TOTAL(5)  + sum_masas*x_cg*z_cg;
I_TOTALCG(6)= I_TOTAL(6)  + sum_masas*y_cg*z_cg;
format longg

OUTPUT=[sum_masas,CG(1),CG(2),CG(3),I_TOTALCG]';
%end