# section.py - CORRECTED TO MATCH MATLAB
# Cross section geometry of design 0 

import numpy as np

geometry = {
    "top_flange_width": 100.0,
    "top_flange_thickness": 1.27,

    "bottom_flange_width": 80.0,
    "bottom_flange_thickness": 1.27,

    "web_thickness": 1.27,
    "web_height": 75.0,  # TOTAL web height (not clear height)

    "side_flange_width": 5.0,    
    "tab_thickness": 1.27,
    
    # Required for buckling calculations:
    "b_in": 80.0 - 2*1.27,  # 77.46 mm - clear width between webs (Case 1)
    "b_out": (100.0 - 80.0) / 2.0  # 10.0 mm - overhang (Case 2)
}

def get_geometry():
    return geometry


def compute_section_properties(geom):
    """
    Compute section properties matching MATLAB logic exactly
    """
    bf_top = geom["top_flange_width"]
    tf_top = geom["top_flange_thickness"]

    bf_bot = geom["bottom_flange_width"]
    tf_bot = geom["bottom_flange_thickness"]

    t_web  = geom["web_thickness"]
    h_web  = geom["web_height"]  # TOTAL height

    tab_w  = geom["side_flange_width"]

    # Calculate CLEAR web height (matching MATLAB)
    clear_web_h = h_web - t_web

    # Y-levels (matching MATLAB logic)
    y_bot_top  = tf_bot                    # top of bottom flange
    y_web_top  = tf_bot + clear_web_h      # top of webs
    y_top_top  = h_web + tf_top            # very top surface

    # Areas + centroids
    # Bottom flange
    A_bot = bf_bot * tf_bot
    y_bot_c = tf_bot / 2.0

    # Side webs (using CLEAR height)
    A_webs = 2 * (t_web * clear_web_h)
    y_webs_c = tf_bot + clear_web_h / 2.0

    # Glue tabs - SAME LEVEL AS TOP FLANGE (matching MATLAB)
    A_tabs = 2 * (tab_w * t_web)
    y_tabs_c = h_web + tf_top / 2.0  # Same as top flange centroid

    # Top flange
    A_top = bf_top * tf_top
    y_top_c = h_web + tf_top / 2.0

    # Neutral axis
    A_total = A_bot + A_webs + A_tabs + A_top

    ybar = (
        A_bot  * y_bot_c
        + A_webs * y_webs_c
        + A_tabs * y_tabs_c
        + A_top  * y_top_c
    ) / A_total

    # Moment of inertia
    I_bot = (bf_bot * tf_bot**3) / 12.0 + A_bot * (y_bot_c - ybar)**2

    I_web_each = (t_web * clear_web_h**3) / 12.0 + (t_web * clear_web_h) * (y_webs_c - ybar)**2
    I_webs = 2 * I_web_each

    I_tab_each = (tab_w * t_web**3) / 12.0 + (tab_w * t_web) * (y_tabs_c - ybar)**2
    I_tabs = 2 * I_tab_each

    I_top = (bf_top * tf_top**3) / 12.0 + A_top * (y_top_c - ybar)**2

    I_total = I_bot + I_webs + I_tabs + I_top

    # Extreme fiber distances
    y_top_dist = y_top_top - ybar
    y_bot_dist = ybar

    # Q at centroid (for shear through webs)
    Q_centroid = 0.0

    # Top flange contribution
    Q_centroid += A_top * (y_top_c - ybar)

    # Tabs contribution (always above NA for this geometry)
    Q_centroid += A_tabs * (y_tabs_c - ybar)

    # Webs portion above NA
    if ybar < y_web_top:
        h_above = y_web_top - ybar
        if h_above > clear_web_h:
            h_above = clear_web_h
        A_web_above = 2 * (t_web * h_above)
        y_web_above_c = ybar + h_above / 2.0
        Q_centroid += A_web_above * (y_web_above_c - ybar)

    # Q for glue line (top flange only)
    Q_glue = A_top * (y_top_c - ybar)

    # Shear widths
    t_shear = 2 * t_web      # two side webs carry shear
    t_glue  = 2 * tab_w      # two glue tabs

    return ybar, y_top_dist, y_bot_dist, I_total, Q_centroid, Q_glue, t_shear, t_glue

