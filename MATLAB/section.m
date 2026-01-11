function sec = section()
% this is the cross-section for design 0
% dfines span, material, geometry, capacities, and section properties.

    %% span and material
    sec.L  = 1200;     % span [mm]
    sec.E  = 4000;     % MPa (N/mm^2), example matboard E
    sec.mu = 0.2;      % Poisson's ratio

    %% Cross-section geometry 
    geom.top_flange_width      = 100;   % mm
    geom.top_flange_thickness  = 1.27;  % mm

    geom.bottom_flange_width   = 80;    % mm
    geom.bottom_flange_thickness = 1.27;

    geom.web_height            = 75;    % mm (overall web depth)
    geom.web_thickness         = 1.27;  % mm

    geom.side_flange_width     = 5;     % mm (glue tab width)

    % train spacing
    sec.spacing_A = 176;   % mm
    sec.spacing_B = 164;   % mm

    % section properties
    sp = compute_section_properties_like_python(geom);

    sec.ybar   = sp.ybar;
    sec.y_top  = sp.y_top;      % ybar -> top surface
    sec.y_bot  = sp.y_bot;      % ybar -> bottom surface
    sec.I      = sp.I;
    sec.Q_web  = sp.Q_centroid;
    sec.Q_glue = sp.Q_glue;
    sec.t_web  = sp.t_shear;
    sec.t_glue = sp.t_glue;

    %% matieral properites (given)
    tension_cap        = 30;   % flexural tension capacity
    compression_cap    = 6;    % flexural compression capacity
    shear_cap_matboard = 4;    % shear in matboard
    shear_cap_glue     = 2;    % shear in glue

    % local flange buckling
    E  = sec.E;
    mu = sec.mu;

    overhang = (geom.top_flange_width - geom.bottom_flange_width) / 2;
    k1 = 0.425;
    sigma_overhang = k1 * pi^2 * E / (12*(1 - mu^2)) * ...
                     (geom.top_flange_thickness / overhang)^2;

    between = geom.bottom_flange_width - 2 * geom.web_thickness;
    k1b = 4.0;
    sigma_between = k1b * pi^2 * E / (12*(1 - mu^2)) * ...
                    (geom.top_flange_thickness / between)^2;

    sigma_flange_crit = min(sigma_overhang, sigma_between);

    % web buckling
    k2 = 6.0;
    sigma_web_crit = k2 * pi^2 * E / (12*(1 - mu^2)) * ...
                     (geom.web_thickness / geom.web_height)^2;

    % shear buckling in web (with diaphragm spacing)
    n_diaphragms = 10;
    spacing_d    = sec.L / (n_diaphragms + 1);

    k_s = 5.0;
    tau_buckling = k_s * pi^2 * E / (12*(1 - mu^2)) * ...
        ((geom.web_thickness / geom.web_height)^2 + ...
         (geom.web_thickness / spacing_d)^2);

    % pack capacities into struct
    cap.tension         = tension_cap;
    cap.compression     = compression_cap;
    cap.flange_buckling = sigma_flange_crit;
    cap.web_buckling    = sigma_web_crit;
    cap.shear_board     = shear_cap_matboard;
    cap.shear_glue      = shear_cap_glue;
    cap.shear_buckling  = tau_buckling;

    sec.cap = cap;
end


% some helper functions
function sp = compute_section_properties_like_python(geom)
    bf_top = geom.top_flange_width;
    tf_top = geom.top_flange_thickness;

    bf_bot = geom.bottom_flange_width;
    tf_bot = geom.bottom_flange_thickness;

    tw     = geom.web_thickness;
    h_web  = geom.web_height;

    side_w = geom.side_flange_width;

    % bottom flange
    A_b = bf_bot * tf_bot;
    y_b = tf_bot / 2;

    % clear web height (between flanges)
    clear_web_h = h_web - tw;

    % two webs
    A_w = 2 * clear_web_h * tw;
    y_w = tf_bot + clear_web_h / 2;

    % two side tabs at top
    A_tabs = 2 * side_w * tw;
    y_tabs = h_web + tf_top / 2;

    % top flange
    A_t = bf_top * tf_top;
    y_t = h_web + tf_top / 2;

    % centroid
    A_total = A_b + A_w + A_tabs + A_t;
    ybar    = (A_b*y_b + A_w*y_w + A_tabs*y_tabs + A_t*y_t) / A_total;

    % second moment of area about NA
    I_b = (bf_bot * tf_bot^3)/12 + A_b * (y_b - ybar)^2;
    I_w = 2 * ((tw * clear_web_h^3)/12 + (clear_web_h * tw) * (y_w - ybar)^2);
    I_tabs = 2 * ((side_w * tw^3)/12 + (side_w * tw) * (y_tabs - ybar)^2);
    I_t = (bf_top * tf_top^3)/12 + A_t * (y_t - ybar)^2;

    I_total = I_b + I_w + I_tabs + I_t;

    % distances to extreme fibres
    y_top_dist = (h_web + tf_top) - ybar;
    y_bot_dist = ybar;

    % Q for web shear (first moment of area above neutral axis)
    % This includes: top flange + side tabs + portion of web above NA
    Q_top_flange = A_t * (y_t - ybar);
    Q_tabs = A_tabs * (y_tabs - ybar);
    
    % Calculate web contribution above NA
    web_top = tf_bot + clear_web_h;
    if ybar < web_top
        % NA is within or below the web
        web_above_h = web_top - ybar;
        if web_above_h > clear_web_h
            web_above_h = clear_web_h;
        end
        % First moment: Area * distance from NA to centroid of that area
        web_cg_above = ybar + web_above_h / 2;
        Q_web_part = 2 * tw * web_above_h * (web_cg_above - ybar);
    else
        % NA is above web (unlikely but handle it)
        Q_web_part = 0;
    end
    
    Q_centroid = Q_top_flange + Q_tabs + Q_web_part;

    % Q for glue line (first moment of top flange about NA)
    % The glue connects top flange to webs
    Q_glue = Q_top_flange;

    % shear thicknesses
    t_shear = 2 * tw;      % two webs carry shear
    t_glue  = 2 * side_w;  % two glue lines (both sides)

    sp = struct( ...
        'ybar',      ybar, ...
        'y_top',     y_top_dist, ...
        'y_bot',     y_bot_dist, ...
        'I',         I_total, ...
        'Q_centroid',Q_centroid, ...
        'Q_glue',    Q_glue, ...
        't_shear',   t_shear, ...
        't_glue',    t_glue);
end