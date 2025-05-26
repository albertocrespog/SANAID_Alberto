function [n] = get_propeller_models(V,Prop_data,T_eng,rho,p_ct,p_ct_red)

% % Cargar los modelos polinómicos
% load('modelos_ct.mat', 'modelos_ct');
% load('modelos_cp.mat', 'modelos_cp');
% load('modelos_ct_red.mat', 'modelos_ct_red');
% load('modelos_cp_red.mat', 'modelos_cp_red');

% Ajuste de modelos polinómicos
%n_curvas_cp = length(modelos_ct);
%grados_pol = 4; % Grado del polinomio
%modelos_cp = cell(1, n_curvas_cp);

% % Ángulos de hélice disponibles
% angulos_helice = [15, 20, 25, 30, 35, 40, 45];
% 
% % Encontrar el índice correspondiente al ángulo de hélice
% indice_curva = find(angulos_helice == beta_pitch);
% if isempty(indice_curva)
%     error('Ángulo de hélice no válido. Debe ser uno de los siguientes: [15, 20, 25, 30, 35, 40, 45]');
% end
% 
% % Obtener el modelo polinómico correspondiente
% p_ct = modelos_ct{indice_curva};
% % p_cp = modelos_cp{indice_curva};
% p_ct_red = modelos_ct_red{indice_curva};
% % p_cp_red = modelos_ct_red{indice_curva};

n = Selection_J_CT_F_des_v2(V,Prop_data,T_eng,rho,p_ct,p_ct_red);
RPM = n*60;

% % Estimacion iniciali apara inicialización
% D = Prop_data.D_prop;
% n_eng = Prop_data.n_eng;
% CT0 = p_ct_red(3);
% CT1 = p_ct_red(2);
% CT2 = p_ct_red(1);
% 
% % Ajuste polinómico 
% N_order_CT = length(p_ct);
% CT_Polyfit = p_ct;
% % Calculates the Advanced Parameter Ratio
% J = @(n) V/(n*D);
% % Determines CT as a funcion of approximation polynomial
% Ct_total =  @(n) CT_Polyfit(N_order_CT+1);
% for j=1:N_order_CT
%     Ct_intermediate = @(n) CT_Polyfit(j)*J(n).^(N_order_CT+1-j);
%     Ct_total = @(n) (Ct_total(n) + Ct_intermediate(n));
% end
% CT = @(n) Ct_total(n);
% 
% % Initialization
% n_sel = -(1/2)*(CT1*V*n_eng*rho*D-sqrt(-4*CT0*CT2*D^2*V^2*n_eng^2*rho^2 + CT1^2*D^2*V^2*n_eng^2*rho^2 + ...
%     4*CT0*T_eng*n_eng*rho))/(CT0*n_eng*rho*D^2);
% 
% % J = @(n) V/(n*D);
% % CT = @(n) CT0 + CT1*J(n) + CT2*J(n)^2;
% Ti_eng = @(n) (n_eng*CT(n)*rho*(n^2)*(D^4));  
% n = fzero(@(n) (T_eng - Ti_eng(n)),n_sel);