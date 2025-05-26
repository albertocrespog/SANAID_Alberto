% Variation of XCG associates to variation of XCG
function x_XCG = Calc_Xcg_vs_W(w_T0,x_XCG,Weight_tier)

% Calculation of xCG for the conditions with/without payload and with
% different batteries

m_TOW = Weight_tier.m_TOW;
m_batteries = Weight_tier.m_batteries;
m_payload = Weight_tier.m_payload;

% Condition 
W(1) = m_TOW;
x_XCG(1) = x_XCG;

% Condition 
W(2) = m_TOW - m_payload/2;
x_XCG(2) = x_XCG;

% Condition 
W(3) = m_TOW - m_payload;
x_XCG(3) = x_XCG;

x_XCG = interp1(W(1:3),Xcg(1:3),m_TOW,'linear');