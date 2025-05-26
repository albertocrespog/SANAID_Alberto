function upw = upwash_calc(AR, X_w, X_c, c_m)
%%  -----------------------------------------------------------------------
%   Academic Stability 2.0
%   28 de Julio de 2015
%   Creado por Álvaro Fernández Cobo
%%   ----------------------------------------------------------------------
%Figura 8.67 Airplane Design, Jan Roskam (PDF 1964)
l_c = abs(X_w-X_c)/c_m;
if l_c > 2
    l_c = 2;
end

% AR_pos = [4 6 9 12];
% [~,i] = min(abs(AR_pos-AR));
% switch i
%     case 1
%         data = load('wing_upw_grad_AR4.dat');
%     case 2
%         data = load('wing_upw_grad_AR6.dat');
%     case 3
%         data = load('wing_upw_grad_AR9.dat');
%     case 4
%         data = load('wing_upw_grad_AR12.dat');
% end
% 
% upw = interp1(data(:,1),data(:,2),l_c,'pchip');

%%  -----------------------------------------------------------------------
%   MODIFIED EMERGENTIA
%   22 de Marzo de 2020
%   Creado por Sergio Esteban
%%   ----------------------------------------------------------------------
%   ----TABLES FOR UPWASH GRADIENT AT PLANE OF SYMMETRY FOR UNSWEPT
%   ----WINGS. FIGURE 4.4.1-73

Y = [4.,6.,9.,12];
X = [.25,.4,.5,.6,.7,.8,1.0,1.2,1.6,2.0];
Z = [1.08,.545,.40,.31,.24,.20,.13,.10,.06,.04; 
      1.18,.680,.52,.40,.32,.27,.19,.15,.10,.08; 
      1.30,.81,.62,.49,.40,.34,.25,.20,.13,.12; 
      1.4,.88,.70,.56,.445,.39,.30,.24,.165,.14];

Yq = AR;
if Yq > 12;
    Yq = 12;
elseif Yq < 4
    Yq = 4;
end
Xq = l_c;
upw = interp2(X,Y,Z,Xq,Yq);

end