function x_XCG = XCG_Estimation_TARSIS(m_TOW,conf_missiles)

Weight_mission = T120_weightconfiguration;

% switch conf_missiles
%     case 0
%     case 1
%     case 2
%     case 3
%     case 4

x(1) = Weight_mission.x_XCG_vec_f(conf_missiles);
x(2) = Weight_mission.x_XCG_vec_f2(conf_missiles);
x(3) = Weight_mission.x_XCG_vec_nf(conf_missiles);

w(1) = Weight_mission.W_vec_f(conf_missiles);
w(2) = Weight_mission.W_vec_f2(conf_missiles);
w(3) = Weight_mission.W_vec_nf(conf_missiles);

x_XCG = interp1(w,x,m_TOW,'spline');

