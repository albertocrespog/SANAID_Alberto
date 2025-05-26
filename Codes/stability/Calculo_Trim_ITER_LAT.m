function [Trim_ITER_LAT,Fig] = Calculo_Trim_ITER_LAT(Geo_calc_ITER,conv_UNITS,V,CL,rho,Stab_Der,CG,m_TO_1,Fig,conditions)

g = conv_UNITS.g;
R2D = conv_UNITS.R2D;
D2R = conv_UNITS.D2R;
         
S_w = Geo_calc_ITER.S_w;
b_w = Geo_calc_ITER.b_w;

q = 0.5*rho*V^2;
         
%beta = 11.5*D2R;
direction_beta = conditions.conditions_TRIM_lat.direction_beta; 
beta = conditions_TRIM_lat.beta;
beta = direction_beta*beta; % corrects direction of beta
beta_vec = conditions_TRIM_lat.beta_vec;
beta_vec = direction_beta*conditions_TRIM_lat.beta_vec; % corrects direction of beta

% gamma = 10*D2R;

T = Stab_Der.T;
CD = Stab_Der.CD;
D = q*S_w*CD;
gamma = asin((T-D)/(m_TO_1*g));
gamma_deg = gamma*R2D;

Trim_ITER_LAT_Viraje.gamma_deg = gamma_deg;

Cyb = Stab_Der.Cyb;
Clb = Stab_Der.Clb;
Cnb = Stab_Der.Cnb;

Cydeltaa = Stab_Der.Cydeltaa;
Cldeltaa = Stab_Der.Cldeltaa;
Cndeltaa = Stab_Der.Cndeltaa;

Cydeltar = Stab_Der.Cydeltar;
Cldeltar = Stab_Der.Cldeltar;
Cndeltar = Stab_Der.Cndeltar;

a1 = Cyb;
a2 = Cydeltaa;
a3 = Cydeltar;
a4 = -m_TO_1*g*cos(gamma)/(q*S_w);

b1 = Clb;
b2 = Cldeltaa;
b3 = Cldeltar;
b4 = 0;

c1 = Cnb;
c2 = Cndeltaa;
c3 = Cndeltar;
c4 = 0;

beta_cnst = beta;
beta = beta_cnst;
deltaa = -(c3*b1*beta+c4*b3-b3*c1*beta-c3*b4)/...
    (-b3*c2+b2*c3);
deltar = (b2*c4-b2*c1*beta-b4*c2+c2*b1*beta)/...
    (-b3*c2+b2*c3);
phi = asin((b2*c3*a1*beta+b2*c4*a3-b2*a3*c1*beta-b4*c2*a3-c2*b3*a1*beta+b4*a2*c3+...
    a2*b3*c1*beta+a3*c2*b1*beta-c3*a2*b1*beta-b3*a2*c4)/(a4*(-b3*c2+b2*c3)));

deltaa_deg = deltaa*R2D;
deltar_deg = deltar*R2D;
phi_deg = phi*R2D;

% Trim_ITER_LAT.deltaa = deltaa;
% Trim_ITER_LAT.deltar = deltar;
% Trim_ITER_LAT.phi = phi;

Trim_ITER_LAT.deltaa_deg = deltaa_deg;
Trim_ITER_LAT.deltar_deg = deltar_deg;
Trim_ITER_LAT.phi_deg = phi_deg;
    
% Variable Beta study
% beta_vec = linspace(0.1*D2R,25*D2R,20);
for i=1:length(beta_vec)
    beta = beta_vec(i);
    deltaa_var(i) = -(c3*b1*beta+c4*b3-b3*c1*beta-c3*b4)/...
        (-b3*c2+b2*c3);
    deltar_var(i) = (b2*c4-b2*c1*beta-b4*c2+c2*b1*beta)/...
        (-b3*c2+b2*c3);
    phi_var(i) = asin((b2*c3*a1*beta+b2*c4*a3-b2*a3*c1*beta-b4*c2*a3-c2*b3*a1*beta+b4*a2*c3+...
        a2*b3*c1*beta+a3*c2*b1*beta-c3*a2*b1*beta-b3*a2*c4)/(a4*(-b3*c2+b2*c3)));
    
    deltaa_deg_var(i) = deltaa_var(i)*R2D;
    deltar_deg_var(i) = deltar_var(i)*R2D;
    phi_deg_var(i) = phi_var(i)*R2D;
end

% Trim_ITER_LAT.deltaa_var = deltaa_var;
% Trim_ITER_LAT.deltar_var = deltar_var;
% Trim_ITER_LAT.phi_var = phi_var;
% 
% Trim_ITER_LAT.deltaa_deg_var = deltaa_deg_var;
% Trim_ITER_LAT.deltar_deg_var = deltar_deg_var;
% Trim_ITER_LAT.phi_deg_var = phi_deg_var;

FIGS=1;
SAVE_FIGS=0;
LS = 2;
TS = 10;

if FIGS ==1
    Fig = Fig + 1;
    figure(Fig)
    plot(beta_vec*R2D,deltaa_deg_var,'b','LineWidth', LS)
    hold on
    plot(beta_vec*R2D,deltar_deg_var,'r','LineWidth', LS)
    plot(beta_vec*R2D,phi_deg_var,'g','LineWidth', LS)
    plot(beta_cnst*R2D,deltaa_deg,'bo','LineWidth', LS*1.5)
    plot(beta_cnst*R2D,deltar_deg,'ro','LineWidth', LS*1.5)
    plot(beta_cnst*R2D,phi_deg,'go','LineWidth', LS*1.5)
    hold off
    title('Lateral Trim Conditions - Sidewind')
    xlabel('\beta (deg)')
    ylabel('trims angles (deg)')
    h_legend=legend('\delta_a','\delta_r','\phi','\delta_{a-TO}','\delta_{r-TO}','\phi_{TO}');
    set(h_legend, 'Location','Best','FontSize',TS)
    grid on
    if SAVE_FIGS==1
        prefix = strcat('Cefiro3_Trim_lateral');
        name   = strcat(prefix,'da_dr_phi');
        saveas(gcf,name,'fig');
        saveas(gca,name,'epsc');
        saveas(gcf,name,'pdf');
        saveas(gcf,name,'bmp');
        saveas(gcf,name,'png');
    end
end