import numpy as np
import matplotlib.pyplot as plt
from math import pi

# ============================================================================
# PARAMETERS
# ============================================================================
L = 1200
spacing_A = 176
spacing_B = 164
TOTAL_LOAD = 400  # total load in Newtons
num_loads = 6
load_per_point = TOTAL_LOAD / num_loads
loads_magnitudes = [load_per_point] * num_loads

slide_step_mm = 1.0
num_points = 1201
x = np.linspace(0, L, num_points)

# Build load positions
offsets = np.array([0, spacing_A, spacing_A + spacing_B, 
                    spacing_A + spacing_B + spacing_A,
                    spacing_A + spacing_B + spacing_A + spacing_B,
                    spacing_A + spacing_B + spacing_A + spacing_B + spacing_A])
total_group_length = offsets[-1]
start_positions = np.arange(0, L - total_group_length + 1e-9, slide_step_mm)

# SFD AND BMD FUNCTIONS
def compute_SFD(load_positions, load_magnitudes):
    total_load = sum(load_magnitudes)
    moment_A = sum(pos * mag for pos, mag in zip(load_positions, load_magnitudes))
    RB = moment_A / L
    RA = total_load - RB
    V = np.ones_like(x) * RA
    for pos, mag in zip(load_positions, load_magnitudes):
        V[x >= pos] -= mag
    return V

def compute_BMD_from_SFD(V):
    dx = x[1] - x[0]
    M = np.cumsum(V) * dx
    M[-1] = 0
    return M

# COMPUTE ENVELOPES
envelope_M_pos = np.full_like(x, -np.inf)
envelope_M_neg = np.full_like(x, np.inf)
envelope_V_pos = np.full_like(x, -np.inf)
envelope_V_neg = np.full_like(x, np.inf)

sample_indices = [0, len(start_positions)//4, len(start_positions)//2, 
                  3*len(start_positions)//4, len(start_positions)-1]
samples_M, samples_V = [], []

print(f"Computing envelopes for {len(start_positions)} load positions...")
for idx, p0 in enumerate(start_positions):
    load_positions = p0 + offsets
    V = compute_SFD(load_positions, loads_magnitudes)
    M = compute_BMD_from_SFD(V)
    
    envelope_M_pos = np.maximum(envelope_M_pos, M)
    envelope_M_neg = np.minimum(envelope_M_neg, M)
    envelope_V_pos = np.maximum(envelope_V_pos, V)
    envelope_V_neg = np.minimum(envelope_V_neg, V)
    
    if idx in sample_indices:
        samples_M.append(M.copy())
        samples_V.append(V.copy())

# Your BME and SFE (unchanged)
envelope_M = envelope_M_pos
envelope_V = envelope_V_pos

# Aliases used in FOS calculation (no change to how envelopes are computed)
BMD_envelope = envelope_M
SFD_envelope = envelope_V

print("Envelope computation complete.")

# ============================================================================
# CROSS-SECTION PROPERTIES (matching your hand calculations)
# ============================================================================
TOP_FLANGE_WIDTH = 100
TOP_FLANGE_THICKNESS = 1.27
SIDE_FLANGE_WIDTH = 5
WEB_HEIGHT = 75
WEB_THICKNESS = 1.27
BOTTOM_FLANGE_WIDTH = 80
BOTTOM_FLANGE_THICKNESS = 1.27

def compute_cross_section_properties():
    """
    Compute cross-section properties matching the exact hand calculation
    """
    t = WEB_THICKNESS
    
    # Component 1: Bottom flange (80mm × 1.27mm)
    A_bottom = BOTTOM_FLANGE_WIDTH * t
    y_bottom = t / 2
    
    # Component 2: Two outer vertical webs (75 - 1.27 = 73.73mm high each)
    web_height = WEB_HEIGHT - t
    A_webs = 2 * web_height * t
    y_webs = t + web_height / 2  # From bottom to center of web
    
    # Component 3: Two side tabs (5mm × 1.27mm each) at top of webs
    A_tabs = 2 * SIDE_FLANGE_WIDTH * t
    y_tabs = WEB_HEIGHT + (TOP_FLANGE_THICKNESS / 2)
    
    # Component 4: Top flange (100mm × 1.27mm)
    A_top = TOP_FLANGE_WIDTH * t
    y_top = WEB_HEIGHT + (TOP_FLANGE_THICKNESS / 2)
    
    # Total area
    total_A = A_bottom + A_webs + A_tabs + A_top
    
    # Calculate ybar (centroid from bottom)
    numerator = (A_bottom * y_bottom + 
                 A_webs * y_webs + 
                 A_tabs * y_tabs + 
                 A_top * y_top)
    ybar = numerator / total_A
    
    # Calculate I using parallel axis theorem
    I_bottom = (BOTTOM_FLANGE_WIDTH * t**3) / 12 + A_bottom * (y_bottom - ybar)**2
    I_webs = 2 * ((t * web_height**3) / 12 + (web_height * t) * (y_webs - ybar)**2)
    I_tabs = 2 * ((SIDE_FLANGE_WIDTH * t**3) / 12 + (SIDE_FLANGE_WIDTH * t) * (y_tabs - ybar)**2)
    I_top = (TOP_FLANGE_WIDTH * t**3) / 12 + A_top * (y_top - ybar)**2
    
    I_total = I_bottom + I_webs + I_tabs + I_top
    
    # Distances from neutral axis
    ytop = (WEB_HEIGHT + t) - ybar  # Distance from NA to top surface
    ybot = ybar  # Distance from NA to bottom surface
    
    # Q at centroid
    A_top_total = A_top + A_tabs
    d_top = y_top - ybar
    
    web_above_height = (WEB_HEIGHT + t) - ybar - t
    A_webs_above = 2 * t * web_above_height
    d_webs = web_above_height / 2
    
    Q_centroid = A_top_total * d_top + A_webs_above * d_webs
    
    # Q at glue joint (top flange)
    A_glue = TOP_FLANGE_WIDTH * t
    d_glue = y_top - ybar
    Q_glue = A_glue * d_glue
    
    # Thicknesses for shear
    t_shear = 2 * t  # Two outer webs
    t_glue = 2 * SIDE_FLANGE_WIDTH  # Two glue joints (side tabs)
    
    return ybar, ytop, ybot, I_total, Q_centroid, Q_glue, t_shear, t_glue

ybar, ytop, ybot, I, Q_centroid, Q_glue, t_shear, t_glue = compute_cross_section_properties()

print(f"\n{'='*60}")
print(f"CROSS-SECTION PROPERTIES")
print(f"{'='*60}")
print(f"ȳ (ybar): {ybar:.3f} mm")
print(f"y_top: {ytop:.3f} mm")
print(f"y_bot: {ybot:.3f} mm")
print(f"I: {I:.6e} mm⁴ = {I/1e6:.6f} × 10⁶ mm⁴")
print(f"Q_centroid: {Q_centroid:.3f} mm³")
print(f"Q_glue: {Q_glue:.3f} mm³")
print(f"t_shear: {t_shear:.2f} mm")
print(f"t_glue: {t_glue:.2f} mm")
print(f"{'='*60}\n")

# ============================================================================
# CALCULATE FOS FOR EACH FAILURE MODE
# ============================================================================
def calculate_FOS_arrays():
    """
    FOS calculation matching iteration1.py:
    - Flexural tension
    - Flexural compression
    - Top flange buckling (overhang & between webs -> min)
    - Web buckling
    - Shear
    - Glue shear
    - Shear buckling
    """
    # Material / plate properties
    E  = 4000   # MPa
    mu = 0.2

    stress_tension_max     = 30  # MPa
    stress_compression_max = 6   # MPa
    tau_matboard_max       = 4   # MPa
    tau_glue_max           = 2   # MPa

    epsilon = 1e-10

    # --- Applied stresses along the bridge (can be + or -) ---

    # Flexural: sigma = M * y / I
    stress_top    = BMD_envelope * ytop / I     # compression at top
    stress_bottom = BMD_envelope * ybot / I     # tension at bottom

    # Shear: tau = V Q / (I b)
    tau_centroid = SFD_envelope * Q_centroid / (I * t_shear)
    tau_glue     = SFD_envelope * Q_glue     / (I * t_glue)

    # Use absolute values for demand in FOS (capacity / |demand|)
    stress_top_abs      = np.abs(stress_top)
    stress_bottom_abs   = np.abs(stress_bottom)
    tau_centroid_abs    = np.abs(tau_centroid)
    tau_glue_abs        = np.abs(tau_glue)

    # --- Buckling capacities (same formulas as iteration1.py) ---

    # Top flange overhang buckling (cantilever plate: one free edge)
    unsupported_flange_overhang = (TOP_FLANGE_WIDTH - BOTTOM_FLANGE_WIDTH) / 2.0
    k1 = 0.425
    sigma_buck_flange_overhang = (
        k1 * pi**2 * E /
        (12 * (1 - mu**2)) *
        (TOP_FLANGE_THICKNESS / unsupported_flange_overhang)**2
    )

    # Top flange buckling between webs (four supported edges)
    flange_between_webs = BOTTOM_FLANGE_WIDTH - 2 * WEB_THICKNESS
    k1b = 4.0
    sigma_buck_flange_between = (
        k1b * pi**2 * E /
        (12 * (1 - mu**2)) *
        (TOP_FLANGE_THICKNESS / flange_between_webs)**2
    )

    # Critical flange buckling stress = min of the two
    sigma_buck_flange_final = min(
        sigma_buck_flange_overhang,
        sigma_buck_flange_between
    )

    # Web buckling in compression
    k2 = 6.0
    sigma_buck_web = (
        k2 * pi**2 * E /
        (12 * (1 - mu**2)) *
        (WEB_THICKNESS / WEB_HEIGHT)**2
    )

    # Shear buckling of webs
    num_diaphragms = 10
    spacing = L / (num_diaphragms + 1)
    k_shear = 5.0
    tau_buck = (
        k_shear * pi**2 * E /
        (12 * (1 - mu**2)) *
        ((WEB_THICKNESS / WEB_HEIGHT)**2 + (WEB_THICKNESS / spacing)**2)
    )

    # --- Factors of safety (capacity / |demand|) ---

    FOS_tension     = stress_tension_max      / (stress_bottom_abs   + epsilon)
    FOS_compression = stress_compression_max / (stress_top_abs       + epsilon)
    FOS_buck_flange = sigma_buck_flange_final / (stress_top_abs      + epsilon)
    FOS_buck_web    = sigma_buck_web          / (stress_top_abs      + epsilon)
    FOS_shear       = tau_matboard_max       / (tau_centroid_abs     + epsilon)
    FOS_glue        = tau_glue_max           / (tau_glue_abs         + epsilon)
    FOS_shear_buck  = tau_buck               / (tau_centroid_abs     + epsilon)

    FOS_min = np.minimum.reduce([
        FOS_tension,
        FOS_compression,
        FOS_buck_flange,
        FOS_buck_web,
        FOS_shear,
        FOS_glue,
        FOS_shear_buck
    ])

    return {
        'tension':      FOS_tension,
        'compression':  FOS_compression,
        'buck_flange':  FOS_buck_flange,
        'buck_web':     FOS_buck_web,
        'shear':        FOS_shear,
        'glue':         FOS_glue,
        'shear_buck':   FOS_shear_buck,
        'min':          FOS_min,
        'stress_top':   stress_top,
        'stress_bottom':stress_bottom,
        'tau_centroid': tau_centroid,
        'tau_glue':     tau_glue,
        'sigma_buck_flange_final': sigma_buck_flange_final,
        'sigma_buck_web':          sigma_buck_web,
        'tau_buck':                tau_buck
    }


FOS_all = calculate_FOS_arrays()

# ============================================================================
# FIND MINIMUM FOS
# ============================================================================
min_FOS = np.min(FOS_all['min'])
fail_idx = np.argmin(FOS_all['min'])
failure_location = x[fail_idx]
P_failure = TOTAL_LOAD * min_FOS

# Match the modes to the keys you actually have
modes = [
    'Flexural Tension',
    'Flexural Compression',
    'Buckling (Top Flange)',
    'Buckling (Web)',
    'Shear',
    'Glue',
    'Shear Buckling'
]

fos_at_fail = [
    FOS_all['tension'][fail_idx],
    FOS_all['compression'][fail_idx],
    FOS_all['buck_flange'][fail_idx],
    FOS_all['buck_web'][fail_idx],
    FOS_all['shear'][fail_idx],
    FOS_all['glue'][fail_idx],
    FOS_all['shear_buck'][fail_idx]
]

failure_mode = modes[np.argmin(fos_at_fail)]

print(f"\n{'='*60}")
print(f"FAILURE ANALYSIS RESULTS")
print(f"{'='*60}")
print(f"Minimum FOS: {min_FOS:.3f}")
print(f"Failure location: {failure_location:.1f} mm")
print(f"Failure mode: {failure_mode}")
print(f"Applied load: {TOTAL_LOAD:.1f} N")
print(f"Predicted failure load: {P_failure:.1f} N")
print(f"\nFOS at failure location:")
for mode, fos in zip(modes, fos_at_fail):
    print(f"  {mode:.<30} {fos:.3f}")

print(f"\nCritical stresses at failure location:")
print(f"  Applied σ_top: {FOS_all['stress_top'][fail_idx]:.3f} MPa")
print(f"  Applied σ_bot: {FOS_all['stress_bottom'][fail_idx]:.3f} MPa")
print(f"  Applied τ_cent: {FOS_all['tau_centroid'][fail_idx]:.3f} MPa")
print(f"  Applied τ_glue: {FOS_all['tau_glue'][fail_idx]:.3f} MPa")

print(f"\nBuckling critical stresses:")
print(f"  σ_cr (top flange, critical): {FOS_all['sigma_buck_flange_final']:.2f} MPa")
print(f"  σ_cr (web):                  {FOS_all['sigma_buck_web']:.2f} MPa")
print(f"  τ_cr (shear buckling):       {FOS_all['tau_buck']:.2f} MPa")
print(f"{'='*60}\n")

# ============================================================================
# VISUALIZATION
# ============================================================================
plt.style.use('default')

# LARGER figure for better visibility
fig = plt.figure(figsize=(24, 14))
fig.patch.set_facecolor('#f7f9fc')

# Better spacing
gs = fig.add_gridspec(2, 3, hspace=0.32, wspace=0.22,
                      left=0.05, right=0.98, top=0.92, bottom=0.06)

colors = ['#6C5CE7', '#00CEC9', '#0984E3', '#E17055', '#00B894']
env_color_m = '#2D3436'
env_color_v = '#D63031'
fos_color_min = '#27AE60'

from matplotlib.lines import Line2D

def style_axis(ax):
    ax.grid(True, alpha=0.25, linestyle='--', linewidth=0.9)
    ax.set_axisbelow(True)
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    ax.spines['left'].set_linewidth(1.3)
    ax.spines['bottom'].set_linewidth(1.3)

# ------------------------- Plot 1: Bending Moment Envelope -------------------
ax1 = fig.add_subplot(gs[0, 0])
for i, M in enumerate(samples_M):
    ax1.plot(x, M, '--', linewidth=1.3, alpha=0.35, color=colors[i % len(colors)])
ax1.plot(x, envelope_M, linewidth=4.0, color=env_color_m, label='$M_{env}(x)$', zorder=10)
ax1.set_xlabel('Position (mm)', fontsize=14, fontweight='600')
ax1.set_ylabel('Bending Moment (N·mm)', fontsize=14, fontweight='600')
ax1.set_title('Bending Moment Envelope', fontsize=17, fontweight='bold', pad=14)
legend_elements = [
    Line2D([0], [0], color=env_color_m, linewidth=4.0, label='Envelope'),
    Line2D([0], [0], color='gray', linewidth=1.3, linestyle='--', alpha=0.6, label='Sample positions')
]
ax1.legend(handles=legend_elements, loc='upper right', framealpha=0.95, fontsize=12)
ax1.set_xlim([0, L])
ax1.tick_params(labelsize=12)
style_axis(ax1)
ax1.axvline(failure_location, linestyle=':', linewidth=2.5, color='#E74C3C', alpha=0.8)

# ------------------------- Plot 2: Shear Force Envelope ----------------------
ax2 = fig.add_subplot(gs[0, 1])
for i, V in enumerate(samples_V):
    ax2.plot(x, V, '--', linewidth=1.3, alpha=0.35, color=colors[i % len(colors)])
ax2.plot(x, envelope_V, linewidth=4.0, color=env_color_v, label='$V_{env}(x)$', zorder=10)
ax2.axhline(y=0, color='#B2BEC3', linestyle='-', linewidth=1.0, alpha=0.7)
ax2.set_xlabel('Position (mm)', fontsize=14, fontweight='600')
ax2.set_ylabel('Shear Force (N)', fontsize=14, fontweight='600')
ax2.set_title('Shear Force Envelope', fontsize=17, fontweight='bold', pad=14)
legend_elements = [
    Line2D([0], [0], color=env_color_v, linewidth=4.0, label='Envelope'),
    Line2D([0], [0], color='gray', linewidth=1.3, linestyle='--', alpha=0.6, label='Sample positions')
]
ax2.legend(handles=legend_elements, loc='upper right', framealpha=0.95, fontsize=12)
ax2.set_xlim([0, L])
ax2.tick_params(labelsize=12)
style_axis(ax2)
ax2.axvline(failure_location, linestyle=':', linewidth=2.5, color='#E74C3C', alpha=0.8)

# ------------------------- Plot 3: Minimum FOS Distribution ------------------
ax3 = fig.add_subplot(gs[0, 2])
ax3.plot(x, FOS_all['min'], linewidth=4.0, color=fos_color_min)
ax3.fill_between(x, 0, FOS_all['min'], color=fos_color_min, alpha=0.12)
ax3.axhline(y=1.0, color='#E74C3C', linestyle='--', linewidth=3.0, label='Failure Threshold', alpha=0.9)
ax3.plot(failure_location, min_FOS, 'o', markersize=15, markeredgewidth=3.0,
         markerfacecolor='#E74C3C', markeredgecolor='#8B0000',
         label=f'Min = {min_FOS:.2f} @ {failure_location:.0f} mm', zorder=15)
ax3.set_xlabel('Position (mm)', fontsize=14, fontweight='600')
ax3.set_ylabel('Factor of Safety', fontsize=14, fontweight='600')
ax3.set_title('Minimum FOS Distribution', fontsize=17, fontweight='bold', pad=14)
ax3.legend(loc='upper right', framealpha=0.95, fontsize=12)
ax3.set_xlim([0, L])
ax3.set_ylim([0, min(10, np.max(FOS_all['min']) * 1.15)])
ax3.tick_params(labelsize=12)
style_axis(ax3)
ax3.axvline(failure_location, linestyle=':', linewidth=2.5, color='#636E72', alpha=0.7)

# ------------------------- Plot 4: FOS Breakdown at Failure ------------------
ax4 = fig.add_subplot(gs[1, 0])
mode_labels = ['Tension', 'Compress.', 'Buck\nFlange', 'Buck\nWeb', 'Shear', 'Glue', 'Shear\nBuck']
bar_colors = ['#27AE60' if v > 2 else '#F39C12' if v > 1 else '#E74C3C' for v in fos_at_fail]
bars = ax4.bar(mode_labels, fos_at_fail, color=bar_colors, alpha=0.9,
               edgecolor='#2D3436', linewidth=1.8, width=0.7)
critical_idx = np.argmin(fos_at_fail)
bars[critical_idx].set_edgecolor('#C0392B')
bars[critical_idx].set_linewidth(4.0)
ax4.axhline(y=1.0, color='#E74C3C', linestyle='--', linewidth=3.0, label='Failure Threshold', alpha=0.9)
ax4.set_ylabel('Factor of Safety', fontsize=14, fontweight='600')
ax4.set_title(f'FOS Breakdown at Failure (x = {failure_location:.0f} mm)', fontsize=17, fontweight='bold', pad=14)
ax4.tick_params(axis='x', rotation=0, labelsize=11)
ax4.tick_params(axis='y', labelsize=12)
ax4.legend(fontsize=12, framealpha=0.95, loc='upper right')
ax4.set_ylim([0, min(15, max(fos_at_fail) * 1.25)])
style_axis(ax4)
for bar, v in zip(bars, fos_at_fail):
    ax4.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.2,
             f'{v:.2f}', ha='center', va='bottom', fontsize=11, fontweight='bold')

# ------------------------- Plot 5: All FOS Curves Distribution ---------------
ax5 = fig.add_subplot(gs[1, 1])

# Define colors for each failure mode (more distinct)
mode_colors = ['#E74C3C', '#3498DB', '#9B59B6', '#1ABC9C', '#F39C12', '#E67E22', '#7F8C8D']

# FOS curves with labels
fos_curves = [
    (FOS_all['tension'], 'Flexural Tension', mode_colors[0]),
    (FOS_all['compression'], 'Flexural Compression', mode_colors[1]),
    (FOS_all['buck_flange'], 'Buckling (Flange)', mode_colors[2]),
    (FOS_all['buck_web'], 'Buckling (Web)', mode_colors[3]),
    (FOS_all['shear'], 'Shear (Matboard)', mode_colors[4]),
    (FOS_all['glue'], 'Shear (Glue)', mode_colors[5]),
    (FOS_all['shear_buck'], 'Shear Buckling', mode_colors[6]),
]

# Plot each FOS curve (clipped to 10 for visibility) - THICKER LINES
for fos_data, label, color in fos_curves:
    fos_clipped = np.clip(fos_data, 0, 10)
    ax5.plot(x, fos_clipped, linewidth=2.8, label=label, color=color, alpha=0.9)

# Failure threshold line
ax5.axhline(y=1.0, color='#C0392B', linestyle='--', linewidth=3.5, label='FOS = 1.0', zorder=10)

# Mark failure location
ax5.axvline(failure_location, linestyle=':', linewidth=2.5, color='#2D3436', alpha=0.7)

ax5.set_xlabel('Position (mm)', fontsize=14, fontweight='600')
ax5.set_ylabel('Factor of Safety', fontsize=14, fontweight='600')
ax5.set_title('All Failure Mode FOS Curves', fontsize=17, fontweight='bold', pad=14)
ax5.set_xlim([0, L])
ax5.set_ylim([0, 10])
ax5.legend(fontsize=10.5, framealpha=0.95, loc='upper right', ncol=2)
ax5.tick_params(labelsize=12)
style_axis(ax5)

# ------------------------- Plot 6: Design Summary ----------------------------
ax6 = fig.add_subplot(gs[1, 2])
ax6.axis('off')

# Status indicator
status = "✗ FAILURE" if min_FOS < 1.0 else "✓ SAFE"
status_color = "#E74C3C" if min_FOS < 1.0 else "#27AE60"

summary_text = f"""
╔══════════════════════════════════════════════╗
║            DESIGN SUMMARY                    ║
╚══════════════════════════════════════════════╝

  Applied Load:        {TOTAL_LOAD} N
  Predicted Failure:   {P_failure:.1f} N
  Factor of Safety:    {min_FOS:.3f}
  Status:              {status}

────────────────────────────────────────────────
  FAILURE INFORMATION
────────────────────────────────────────────────
  Mode:      {failure_mode}
  Location:  {failure_location:.0f} mm

────────────────────────────────────────────────
  CROSS-SECTION PROPERTIES
────────────────────────────────────────────────
  Centroid (ȳ):   {ybar:.2f} mm
  I:              {I/1e6:.4f} × 10⁶ mm⁴

  Top Flange:     {TOP_FLANGE_WIDTH} × {TOP_FLANGE_THICKNESS} mm
  Bottom Flange:  {BOTTOM_FLANGE_WIDTH} × {BOTTOM_FLANGE_THICKNESS} mm
  Web Height:     {WEB_HEIGHT} mm
  Web Thickness:  {WEB_THICKNESS} mm

────────────────────────────────────────────────
  MATERIAL CAPACITIES
────────────────────────────────────────────────
  Tension:         30 MPa
  Compression:      6 MPa
  Matboard Shear:   4 MPa
  Glue Shear:       2 MPa
"""

ax6.text(0.02, 0.98, summary_text, transform=ax6.transAxes,
         fontsize=12, verticalalignment='top', fontfamily='monospace',
         bbox=dict(boxstyle='round,pad=0.8', facecolor='white', alpha=0.95, 
                   edgecolor='#2D3436', linewidth=2.5))

# Large status badge at bottom
ax6.text(0.50, 0.02, status, transform=ax6.transAxes,
         fontsize=22, fontweight='bold', ha='center', va='bottom', color=status_color,
         bbox=dict(boxstyle='round,pad=0.6', facecolor='white', edgecolor=status_color, linewidth=4))

plt.suptitle(f'Bridge Failure Analysis — Predicted Failure: {P_failure:.1f} N at {failure_location:.0f} mm ({failure_mode})',
             fontsize=20, fontweight='bold', y=0.98)

# Save high-res output
plt.savefig('bridge_analysis_output.png', dpi=150, bbox_inches='tight', facecolor='#f7f9fc')
plt.show()

print("\n✓ Plot saved as 'bridge_analysis_output.png'")
