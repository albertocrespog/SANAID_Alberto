function [KVB, KBV] = get_Vtail_body_interference(Body_Geo, x_w2_LE, b_w2_e)

length_x_position = Body_Geo.length_x_position;
width_x_position = Body_Geo.width_x_position;
d = interp1(length_x_position,width_x_position,x_w2_LE);

db= d/b_w2_e;

%Fig. 4.3.1.2-10 DATCOM CAP4-B. 

% C                   FIGURE 4.3.1,2-10 KWB
% C
      d_b = [0.0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0];
      KVeeB = [1.0,1.08,1.16,1.26,1.36,1.46,1.56,1.67,1.78,1.89,2.0];
% C
% C                   FIGURE 4.3.1.2-10 KBW
% C
      KBVee = [0.0,.13,.29,.45,.62,.80,1.0,1.22,1.45,1.70,2.0];
 
KVB = interp1(d_b, KVeeB, db);
KBV = interp1(d_b, KBVee, db);
end
