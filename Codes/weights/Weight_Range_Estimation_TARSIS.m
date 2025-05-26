function w = Weight_Range_Estimation_TARSIS(conf_missiles)

Weight_mission = T120_weightconfiguration;

w(1) = Weight_mission.W_vec_f(conf_missiles);
w(2) = Weight_mission.W_vec_f2(conf_missiles);
w(3) = Weight_mission.W_vec_nf(conf_missiles);
