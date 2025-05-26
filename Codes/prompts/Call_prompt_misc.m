function answer_misc = Call_prompt_misc

prompt = {'Enter Prop Diameter Estimation (in):','Enter scale factor %:','Enter Weight Estimation Method: case 1 - composite Céfiro III, case 2 - wood Céfiro I'};
dlgtitle = 'Input';
dims = [1 35];
opts.Interpreter = 'tex';
SF_in = 1;
SF_str = num2str(SF_in);
% Weight Estimation
Weight_in = 1;
Weight_str = num2str(Weight_in);
% Scaling Factor (in)
D_propin_in = 30;
D_propin_str = num2str(D_propin_in);
definput = {D_propin_str,SF_str,Weight_str};
answer_misc = inputdlg(prompt,dlgtitle,dims,definput);