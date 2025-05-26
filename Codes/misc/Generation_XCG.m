function XCG_DATA = Generation_XCG(Geo_tier,Weight_tier);

Weight_tier.m_TOW = m_TOW;
Weight_tier.m_estructure = m_estructure;
m_w1 = Weight_tier.m_w1;
m_w2 = Weight_tier.m_w2;
m_fus = Weight_tier.m_fus;
m_fairing = Weight_tier.m_fairing;
m_engs = Weight_tier.m_engs;
m_props = Weight_tier.m_props;
m_servos = Weight_tier.m_servos;
m_wiring = Weight_tier.m_wiring;
m_metal = Weight_tier.m_metal;
m_payload = Weight_tier.m_payload;
m_systems = Weight_tier.m_systems;
m_batteries = Weight_tier.m_batteries;

% m_estructure = m_w1 + m_w2 + m_fus + m_fairing + m_engs + m_props + m_servos + m_wiring + m_metal ;
m_estructure = m_w1 + m_w2 + m_fairing + m_engs + m_props + m_servos + m_wiring + m_metal ;
m_payload = 2;
m_systems = 0.5;
m_batteries = 8;

m_TOW = m_estructure + m_payload + m_systems + m_batteries;