import numpy as np

geometry = {
    "top_flange_width": 100.0,
    "top_flange_thickness": 1.27,     # single ply
    "bottom_flange_width": 92.54,
    "bottom_flange_thickness": 1.27,
    "web_thickness": 1.27,
    "web_height": 56.19,              # clear height
    "tab_width": 5.0,
    "n_webs": 4
}

def get_geometry():
    return geometry


def compute_section_properties(geom):

    # pull geometry
    bf_top = geom["top_flange_width"]
    tf_top = geom["top_flange_thickness"]
    bf_bot = geom["bottom_flange_width"]
    tf_bot = geom["bottom_flange_thickness"]
    web_t  = geom["web_thickness"]
    web_h  = geom["web_height"]
    tab_w  = geom["tab_width"]
    n_webs = geom["n_webs"]

    # top flange is double ply
    tf_top_total = 2 * tf_top

 
    y_bot_top  = tf_bot                         # top of bottom flange
    y_web_top  = tf_bot + web_h                 # top of webs = underside of top flange
    y_top_top  = y_web_top + tf_top_total       # very top surface (should equal 60)


    y_bot_c  = tf_bot / 2.0
    y_webs_c = tf_bot + web_h / 2.0


    y_tabs_c = y_web_top - web_t / 2.0


    y_top_c = y_web_top + tf_top_total / 2.0

    A_bot  = bf_bot * tf_bot
    A_webs = n_webs * web_t * web_h
    A_tabs = 2 * (tab_w * web_t)
    A_top  = bf_top * tf_top_total

    A_total = A_bot + A_webs + A_tabs + A_top

    ybar = (
        A_bot  * y_bot_c +
        A_webs * y_webs_c +
        A_tabs * y_tabs_c +
        A_top  * y_top_c
    ) / A_total

    I_bot = (bf_bot * tf_bot**3)/12 + A_bot*(y_bot_c - ybar)**2

    I_web_each = (web_t * web_h**3)/12 + (web_t*web_h)*(y_webs_c - ybar)**2
    I_webs = n_webs * I_web_each

    I_tab_each = (tab_w * web_t**3)/12 + (tab_w*web_t)*(y_tabs_c - ybar)**2
    I_tabs = 2 * I_tab_each

    I_top = (bf_top * tf_top_total**3)/12 + A_top*(y_top_c - ybar)**2

    I_total = I_bot + I_webs + I_tabs + I_top


    y_top_dist = y_top_top - ybar
    y_bot_dist = ybar


    Q_centroid = 0.0

    # top flange
    Q_centroid += A_top * (y_top_c - ybar)

    # tabs
    if ybar < y_web_top:
        Q_centroid += A_tabs * (y_tabs_c - ybar)

    # webs above NA
    if ybar < y_web_top:
        h_above = y_web_top - ybar
        A_webs_above = n_webs * web_t * h_above
        y_webs_above_c = ybar + h_above/2
        Q_centroid += A_webs_above * (y_webs_above_c - ybar)

    # glue-line Q: only top flange above glue line
    Q_glue = A_top * (y_top_c - ybar)

    # shear thicknesses
    t_shear = n_webs * web_t
    t_glue  = 2 * tab_w

    return ybar, y_top_dist, y_bot_dist, I_total, Q_centroid, Q_glue, t_shear, t_glue
