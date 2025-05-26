function Lambda_nth_w = Get_Nth_Lambda(Nth,AR_w,Lambda_LE_w,lambda_w)

% Calculates the Nth sweepback of the nth-chord line as a function of
% geometry
Lambda_nth_w = atan((1/AR_w)*(AR_w*tan(Lambda_LE_w) - 4*(Nth)*((1-lambda_w)/(1+lambda_w))));
