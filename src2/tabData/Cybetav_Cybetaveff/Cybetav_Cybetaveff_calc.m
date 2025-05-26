%%  -----------------------------------------------------------------------
%   EMERGENTIA
%   11 de Diciembre de 2019
%   Creado por Sergio Esteban
%%   ----------------------------------------------------------------------
% Calculo de alphaCL_alphaCl: The wing-fuselage-"equivalent horizontal tail" interference on the airplane sideforce-coefficient-due-to-sideslip
% derivative of equivalent twin vertical tails is obtained from Figure 10.17 in Airplane Design Part VI and is a function of the equivalent 
% twin vertical tail span, the fuselage depth at the quarter chord point of the V-tail, the distance between the two equivalent vertical tails, 
% and the fuselage length:
function Cybeta_v_Cybeta_eff = Cybetav_Cybetaveff_calc(b_h,l_f,r_1,b_v)

Yq = b_h/l_f;
Xq = 2*r_1/b_v;

if Xq > 1
    Xq = 1;
elseif Xq < 0 
        Xq = 0;
end

if Yq > 1
    Yq = 1;
elseif Yq < 0.2 
        Yq = 0.2;
end


% FIGURE 5.3.1.1-24C
Y = [.2 .4 .6 .8 1.0];
X = [0.0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1.0];
Z =[ 1.0,.90,.818,.741,.68,.631,.593,.563,.542,.527,.510;
     1.0,.93,.87,.825,.784,.752,.730,.710,.696,.688,.680;
     1.0,.95,.911,.882,.859,.844,.830,.820,.812,.805,.800;
     1.0,.968,.941,.923,.908,.898,.890,.883,.881,.880,.879;
     1.0,.979,.962,.954,.949,.943,.938,.935,.935,.931,.931];
Cybeta_v_Cybeta_eff = interp2(X,Y,Z,Xq,Yq);
     