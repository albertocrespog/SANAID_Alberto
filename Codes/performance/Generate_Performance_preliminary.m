function Performance_preliminar = Generate_Performance_preliminary(OUTPUT_read_XLSX);

Performance_preliminar = OUTPUT_read_XLSX.Performance_pre_flags;

% Atmospheric conditions
[Temp_init,rho_init,p_init,a_init]=atmos_inter_mio(Performance_preliminar.h);
Performance_preliminar.Temp = Temp_init;
Performance_preliminar.rho = rho_init;
Performance_preliminar.p = p_init;
Performance_preliminar.a = a_init;
Mach_init = Performance_preliminar.V/a_init;
Performance_preliminar.Mach = Mach_init;
q_inf_init = 0.5*rho_init*(Performance_preliminar.V)^2;
Performance_preliminar.q_inf = q_inf_init;
Performance_preliminar.Flight_cruise = OUTPUT_read_XLSX.Performance_pre_flags.Flight_cruise; %
Performance_preliminar.Flight_takeoff = OUTPUT_read_XLSX.Performance_pre_flags.Flight_takeoff; %
Performance_preliminar.Flight_climb = OUTPUT_read_XLSX.Performance_pre_flags.Flight_climb; %
