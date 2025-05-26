function XCG_data = Generation_XCG_data(Geo_input_tier,Prop_data,conv_UNITS,AC_CONFIGURATION,SF)
%% Location of the Xcg
% Defines the XCG location
x_XCG = SF*1.1747; 
y_XCG = SF*0;
z_XCG = SF*0;

% Storing DATA
XCG_data.x_XCG = x_XCG;
XCG_data.y_XCG = y_XCG;
XCG_data.z_XCG = z_XCG;