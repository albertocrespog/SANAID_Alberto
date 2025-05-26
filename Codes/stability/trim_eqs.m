function F = trim_eqs(x, CL_needed, CL0_ac, CL_alpha_ac, CL_delta, CM0_needed, CM_alpha_ac, CM_delta);
F = [ CL_needed - CL0_ac - CL_alpha_ac*x(1) - CL_delta*x(2)
  CM0_needed + CM_alpha_ac*x(1) + CM_delta*x(2)];
end