function sec = section2()
% our final cross section
%  L = 1269 mm, H_total = 60 mm
% our cross section used webs total (each 1.27 mm thick)
% bottom flange is single layer
% Top flange is double layered

    %% span + material
    sec.L  = 1269;      % mm
    sec.E  = 4000;      % MPa
    sec.mu = 0.2;

    %% train spacings
    sec.spacing_A = 176;   % mm
    sec.spacing_B = 164;   % mm

    %% geometry
    g.t        = 1.27;       % base matboard thickness [mm]
    g.b_top    = 100.0;      % top flange width [mm]
    g.b_bot    = 92.54;      % bottom flange width [mm]
    g.tab_w    = 5.0;        % glue tab width each side [mm]
    g.H_total  = 60.0;       % total height [mm]
    g.n_webs   = 4;          % number of webs (each thickness t)

    % flange thicknesses
    g.t_bot = g.t;           % bottom flange = 1 layer
    g.t_top = 2*g.t;         % top flange = 2 layers stacked

    % clear web height between flanges
    g.h_web = g.H_total - g.t_bot - g.t_top;  % mm (clear web height)

    % diaphragm spacing a for Case 4
    g.a_diaphragm = 40.0;    % mm

    % paired webs act as one plate (no spacing)
    g.web_pair_thickness = 2*g.t;  % 2.54 mm

    %% Section properties (MATCHING PYTHON)
    sp = compute_section_properties_python_match(g);

    sec.ybar   = sp.ybar;
    sec.y_top  = sp.y_top;
    sec.y_bot  = sp.y_bot;
    sec.I      = sp.I;
    sec.Q_web  = sp.Q_web;
    sec.Q_glue = sp.Q_glue;

    % Effective shear widths b in tau = VQ/(I b)
    sec.t_web  = sp.b_web_eff;   % total web shear width
    sec.t_glue = sp.b_glue_eff;  % total glue width

    %% Material capacities (MPa)
    tension_cap        = 30;
    compression_cap    = 6;
    shear_cap_matboard = 4;
    shear_cap_glue     = 2;

    E  = sec.E;
    mu = sec.mu;

    % buckling parameters 
    t_flange = g.t_top;                 % thickness of top flange plate (total)
    b_web    = g.web_pair_thickness;    % paired web plate thickness
    h_web    = g.h_web;                 % clear web height
    ytop     = sec.y_top;               % distance NA -> top fibre

    % buckling cases

    % case 1
    b_in = g.b_bot - 2*b_web;  
    k1 = 4.0;
    sigma_case1 = k1*pi^2*E/(12*(1-mu^2)) * (t_flange/b_in)^2;

    % case 2
    b_out = (g.b_top - g.b_bot)/2;
    k2 = 0.4254;
    sigma_case2 = k2*pi^2*E/(12*(1-mu^2)) * (t_flange/b_out)^2;

    sigma_flange_crit = min(sigma_case1, sigma_case2);

    % case 3
    k3 = 6.0;
    sigma_web_crit = k3*pi^2*E/(12*(1-mu^2)) * (b_web/ytop)^2;

    % case 4
    a  = g.a_diaphragm;
    k4 = 5.0;
    tau_buckling = k4*pi^2*E/(12*(1-mu^2)) * ( (b_web/h_web)^2 + (b_web/a)^2 );

    %% material properites + caps
    cap.tension         = tension_cap;
    cap.compression     = compression_cap;
    cap.flange_buckling = sigma_flange_crit;  % governing case 1/2
    cap.web_buckling    = sigma_web_crit;     % case 3
    cap.shear_board     = shear_cap_matboard;
    cap.shear_glue      = shear_cap_glue;
    cap.shear_buckling  = tau_buckling;       % case 4

    sec.cap = cap;
end


% Section properties 
function sp = compute_section_properties_python_match(g)
    t      = g.t;
    bt     = g.b_top;
    tt     = g.t_top;      % top flange total thickness (2t)
    bb     = g.b_bot;
    tb     = g.t_bot;      % bottom flange thickness (t)
    hw     = g.h_web;      % clear web height
    tab_w  = g.tab_w;
    n_webs = g.n_webs;

    % Y-LEVELS (matching Python)
    y_bot_top  = tb;                    % top of bottom flange
    y_web_top  = tb + hw;               % top of webs = underside of top flange
    y_top_top  = y_web_top + tt;        % very top surface (should = 60)

    % Centroids of each component
    y_bot_c    = tb / 2.0;
    y_webs_c   = tb + hw / 2.0;
    
    % Tabs live INSIDE at top, under the double top flange
    y_tabs_c   = y_web_top - t / 2.0;
    
    y_top_c    = y_web_top + tt / 2.0;

    % areas (tabs ARE structural)
    A_bot  = bb * tb;
    A_webs = n_webs * t * hw;
    A_tabs = 2 * (tab_w * t);
    A_top  = bt * tt;

    A_total = A_bot + A_webs + A_tabs + A_top;

    % ybar
    ybar = (A_bot*y_bot_c + A_webs*y_webs_c + A_tabs*y_tabs_c + A_top*y_top_c) / A_total;

    % intertia
    I_bot = (bb*tb^3)/12 + A_bot*(y_bot_c - ybar)^2;
    
    I_web_each = (t*hw^3)/12 + (t*hw)*(y_webs_c - ybar)^2;
    I_webs = n_webs * I_web_each;
    
    I_tab_each = (tab_w*t^3)/12 + (tab_w*t)*(y_tabs_c - ybar)^2;
    I_tabs = 2 * I_tab_each;
    
    I_top = (bt*tt^3)/12 + A_top*(y_top_c - ybar)^2;

    I_total = I_bot + I_webs + I_tabs + I_top;


    % shear free areas (Q)
    y_top_dist = y_top_top - ybar;
    y_bot_dist = ybar;

    % Q at centroid
    Q_web = 0;
    
    Q_web = Q_web + A_top*(y_top_c - ybar);
    
    if ybar < y_web_top
        Q_web = Q_web + A_tabs*(y_tabs_c - ybar);
    end
    
    if ybar < y_web_top
        h_above = y_web_top - ybar;
        A_webs_above = n_webs * t * h_above;
        y_webs_above_c = ybar + h_above/2;
        Q_web = Q_web + A_webs_above*(y_webs_above_c - ybar);
    end

    % Q for glue
    Q_glue = A_top*(y_top_c - ybar);

    b_web_eff  = n_webs * t;     
    b_glue_eff = 2 * tab_w;      



    sp.ybar      = ybar;
    sp.y_top     = y_top_dist;
    sp.y_bot     = y_bot_dist;
    sp.I         = I_total;
    sp.Q_web     = Q_web;
    sp.Q_glue    = Q_glue;
    sp.b_web_eff = b_web_eff;
    sp.b_glue_eff= b_glue_eff;
end
