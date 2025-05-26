function [Stab_Der_parts Stab_Der] = get_tail_incidence_derivatives(AC_CONFIGURATION,modelo,TRIM_RESULTS,Stab_Der,Stab_Der_parts,Trim_ITER,OUTPUT_read_XLSX)

W1 = AC_CONFIGURATION.W1;
VTP = AC_CONFIGURATION.VTP;
HTP = AC_CONFIGURATION.HTP;
Can = AC_CONFIGURATION.Can;
Vee = AC_CONFIGURATION.Vee;
Vee2 = AC_CONFIGURATION.Vee2;
Nac = AC_CONFIGURATION.Nac;

prop_wash_effect = OUTPUT_read_XLSX.Stability_flags.prop_wash_effect;

Stab_Der.CD_Vtail_i = 0;
Stab_Der.CL_Vtail_i = 0;
Stab_Der.CM_Vtail_i = 0;
Stab_Der.Cy_Vtail_i = 0;
Stab_Der.Cl_Vtail_i = 0;
Stab_Der.Cn_Vtail_i = 0;
Stab_Der.CD_HTP_i = 0;
Stab_Der.CL_HTP_i = 0;
Stab_Der.CM_HTP_i = 0;

Stab_Der_parts.CD_Vtail_i = 0;
Stab_Der_parts.CL_Vtail_i = 0;
Stab_Der_parts.CM_Vtail_i = 0;
Stab_Der_parts.Cy_Vtail_i = 0;
Stab_Der_parts.Cl_Vtail_i = 0;
Stab_Der_parts.Cn_Vtail_i = 0;
Stab_Der_parts.CD_HTP_i = 0;
Stab_Der_parts.CL_HTP_i = 0;
Stab_Der_parts.CM_HTP_i = 0;
