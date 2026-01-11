% includes anaylsis similar to the one in python(but not as deep) %
    includes a 3D visual of the bridge as it deflects using MAT2
    % some boundaries were set to avoid drifting

    % note
    : my way of doing the train load probably
      wasn't the smartest but it % works , much better methods can be employed % this might be
      quite heavy to run,
    if you're experiencing issues go to the % variable nsteps and change it to 180 / 120 or whatever you need

    clear;
clc;
close all;

% defining section, a new file was made for each different cross-section
% to make it easy to switch from one cross section to another
SEC_FUNC = @section;
sec = SEC_FUNC();

L = sec.L;
E = sec.E;
I = sec.I;
ybar = sec.ybar;
y_top = sec.y_top;
y_bot = sec.y_bot;
Q_web = sec.Q_web;
Q_glue = sec.Q_glue;
t_web = sec.t_web;
t_glue = sec.t_glue;
caps = sec.cap;

fprintf('\n========== BRIDGE SIMULATION STARTING ==========\n');
fprintf('Span: %.0f mm | E: %.0f MPa | I: %.1f mm^4\n', L, E, I);
fprintf('Total train load: 400 N (6 axles × %.1f N)\n\n', 400 / 6);

n = 600;
x = linspace(0, L, n + 1);

%% Train Load for Loadcase1
total_load = 400;
num_axles = 6;
P_axle = total_load / num_axles;
for
  A = sec.spacing_A;
% these are the two different spacings between the wheels B = sec.spacing_B;

offsets = [ 0, A, A + B, A + B + A, A + B + A + B, A + B + A + B + A ];
train_length = offsets(end);

% % setting up visulization of the beam W_bridge = 100;
H_bridge = 60;
zmin = -W_bridge / 2;
zmax = W_bridge / 2;

y_top_face0 = 0;

% Animation setup nSteps = 260;
front_start = -train_length;
front_end = L + 60;
front_pos = linspace(front_start, front_end, nSteps);

hFig = figure('Name', 'Bridge Structural Analysis', 'NumberTitle', 'off');
clf(hFig);
set(hFig, 'Color', [0.97 0.97 0.97], 'Units', 'normalized', 'Position',
    [0.03 0.08 0.94 0.84]);

ax = subplot('Position', [0.05 0.08 0.58 0.88]);
hold(ax, 'on');
grid(ax, 'on');
box(ax, 'on');
view(ax, 130, 25);
camproj(ax, 'perspective');
xlabel(ax, 'Longitudinal Position (mm)', 'FontSize', 12, 'FontWeight', 'normal',
       'FontName', 'Segoe UI');
ylabel(ax, 'Vertical Deflection (mm)', 'FontSize', 12, 'FontWeight', 'normal',
       'FontName', 'Segoe UI');
zlabel(ax, 'Lateral Position (mm)', 'FontSize', 12, 'FontWeight', 'normal',
       'FontName', 'Segoe UI');
set(ax, 'FontName', 'Segoe UI', 'FontSize', 10.5, 'Color', [0.99 0.99 0.99]);
set(ax, 'GridAlpha', 0.12, 'GridLineStyle', '-', 'LineWidth', 0.6);

% Side panel for anayalsis
ax_text = subplot('Position', [0.66 0.08 0.32 0.88]);
axis(ax_text, 'off');

failed = false;
fail_x = NaN;
fail_mode = '';
FOS_at_fail = NaN;
P_on_fail = NaN;
eps_small = 1e-10;

min_FOS_history = inf;
critical_mode_history = '';
max_deflection_history = 0;

pause_done = false;
% so we only pause once at first failure

        modes = {... 'Flexural Tension',     ... 'Flexural Compression',
                 ... 'Top Flange Buckling',  ... 'Web Buckling',
                 ... 'Shear in Web',         ... 'Shear in Glue',
                 ... 'Shear Buckling in Web'};

nModes = numel(modes);
min_FOS_mode_history = inf(nModes, 1);

% panel layout constants(kept tight so everything fits) panel_left = 0.06;
panel_right = 0.94;
panel_top = 0.955;
panel_bot = 0.045;

for
  k = 1 : nSteps

              front_x = front_pos(k);
axle_x = front_x + offsets;

on_span = (axle_x >= 0) & (axle_x <= L);
load_x = axle_x(on_span);
load_P = P_axle * ones(1, numel(load_x));
P_on = sum(load_P);

if isempty (load_x)
  V = zeros(size(x));
M = zeros(size(x));
RA = 0;
RB = 0;
else[V, M, RA, RB] = compute_SFD_BMD_local(x, L, load_x, load_P);
end

    % defaults so code never breaks when no load FOS_global = inf;
critical_mode = 'None';
fail_location = NaN;

% == == == == == == == == == == == == == == == == == == == == == == == == == ==
    == == ==
    = % FOS CALCULATION(RUNS EVEN AFTER FAILURE) % == == == == == == == == == ==
      == == == == == == == == == == == == == == == == == == ==
    = if P_on > 0 stress_top = -M.*y_top./ I;
stress_bot = M.*y_bot./ I;

stress_bot_tens = max(stress_bot, 0);
% tension only at bottom stress_top_comp = max(-stress_top, 0);
% compression only at top

        tau_web = abs(V).*Q_web./ (I * t_web);
tau_glue = abs(V).*Q_glue./ (I * t_glue);

FOS_tension = caps.tension./ (stress_bot_tens + eps_small);
FOS_compression = caps.compression./ (stress_top_comp + eps_small);
FOS_flange_buck = caps.flange_buckling./ (stress_top_comp + eps_small);
FOS_web_buck = caps.web_buckling./ (stress_top_comp + eps_small);
FOS_shear = caps.shear_board./ (tau_web + eps_small);
FOS_glue = caps.shear_glue./ (tau_glue + eps_small);
FOS_shear_buck = caps.shear_buckling./ (tau_web + eps_small);

FOS_frame_mins = [... min(FOS_tension);... min(FOS_compression);
                  ... min(FOS_flange_buck);... min(FOS_web_buck);
                  ... min(FOS_shear);... min(FOS_glue);... min(FOS_shear_buck)];

% keep updating historical mins no matter what min_FOS_mode_history =
    min(min_FOS_mode_history, FOS_frame_mins);

FOS_stack = cat(1, FOS_tension, FOS_compression, FOS_flange_buck,
                ... FOS_web_buck, FOS_shear, FOS_glue, FOS_shear_buck);

[ FOS_min, mode_idx ] = min(FOS_stack, [], 1);
[ FOS_global, iFail ] = min(FOS_min);
fail_location = x(iFail);
critical_mode = modes{mode_idx(iFail)};

if FOS_global
  < min_FOS_history min_FOS_history = FOS_global;
critical_mode_history = critical_mode;
end

    % only *set *failure once,
    but keep sim running after if ~failed &&FOS_global <= 1.0 failed = true;
fail_x = x(iFail);
fail_mode = critical_mode;
FOS_at_fail = FOS_global;
P_on_fail = P_on;

fprintf('\n!!! BRIDGE FAILURE DETECTED !!!\n');
fprintf('Frame: %d/%d\n', k, nSteps);
fprintf('Location: x = %.1f mm\n', fail_x);
fprintf('Failure mode: %s\n', fail_mode);
fprintf('Load on span: %.1f N (%.1f%% of total)\n', ... P_on_fail,
        100 * P_on_fail / total_load);
fprintf('FOS at failure: %.4f\n\n', FOS_at_fail);
end end

    % Deflection calculation(MAT2) w = zeros(size(x));

if P_on
  > 0 && max(abs(M)) > 1e-10 kappa = M./ (E * I);
p_fit = polyfit(x, kappa, 1);
kappa = kappa - polyval(p_fit, x);

theta = cumtrapz(x, kappa);
theta = theta - (theta(end) / L) * x;

w = cumtrapz(x, theta);
w = w - w(1) - ((w(end) - w(1)) / L) * x;

w(1) = 0;
w(end) = 0;

if max (abs(w))
  > L / 50 w = w * (L / 50) / max(abs(w));
end end

    max_deflection = max(abs(w));

if max_deflection
  > max_deflection_history max_deflection_history = max_deflection;
end

    if max_deflection > 1e-6 if max_deflection >
    L / 100 scale_factor = 0.15 * L / max_deflection;
else scale_factor = 0.25 * L / max_deflection;
end else scale_factor = 1;
end

    w_smooth = w * scale_factor;

% 3D Animation cla(ax);
hold(ax, 'on');
grid(ax, 'on');

sup_h = 50;
draw_support_block(ax, 0, sup_h, W_bridge);
draw_support_block(ax, L, sup_h, W_bridge);

if
  ~failed draw_box_beam_smooth(ax, x, w_smooth, y_top_face0, H_bridge, zmin,
                               zmax, ...[0.90 0.95 0.99], [0.28 0.48 0.78]);
else
  draw_box_beam_smooth(ax, x, w_smooth, y_top_face0, H_bridge, zmin, zmax,
                       ...[1.0 0.72 0.72], [0.87 0.17 0.17]);

if fail_x
  >= min(x) && fail_x <= max(x) y_fail =
      interp1(x, w_smooth, fail_x, 'linear', 0) + y_top_face0;
plot3(ax, fail_x, y_fail, 0, 'rx', 'MarkerSize', 24, 'LineWidth', 5);
text(fail_x, y_fail + 75, 0, '⚠ FAILURE', ... 'FontSize', 13, 'FontWeight',
     'bold', 'Color', [0.92 0 0], ... 'HorizontalAlignment', 'center', 'Parent',
     ax, 'FontName', 'Segoe UI');
        end
    end

    if ~isempty(load_x)
        for j = 1:numel(load_x)
            px = load_x(j);

        if px
          >= min(x) && px <= max(x) y_top_here =
              interp1(x, w_smooth, px, 'linear', 0) + y_top_face0;
        else
          y_top_here = y_top_face0;
        end

            arrow_len = 60;
        quiver3(ax, px, y_top_here + H_bridge + 15, 0, 0, -arrow_len, 0,
                ... 'Color', [0.94 0.14 0.14], 'LineWidth', 3.4, 'MaxHeadSize',
                1.6);
        end end

            xlim(ax, [ -150, L + 150 ]);
        ylim(ax, [ -0.4 * L, 0.3 * L ]);
        zlim(ax, [ -W_bridge * 0.7, W_bridge * 0.7 ]);
        view(ax, 130, 25);

        axis(ax, 'vis3d');
        daspect(ax, [1 1 1]);
        camlight(ax, 'headlight');
        lighting(ax, 'gouraud');

        if
          ~failed title(ax,
                        sprintf('Bridge Structural Model  |  Frame %d/%d', k,
                                nSteps),
                        ... 'FontWeight', 'normal', 'FontSize', 14, 'FontName',
                        'Segoe UI', 'Color', [0.15 0.15 0.15]);
        else
          title(ax,
                sprintf('FAILURE DETECTED  |  x = %.0f mm  |  %s', fail_x,
                        fail_mode),
                ... 'FontWeight', 'bold', 'Color', [0.82 0 0], 'FontSize', 14,
                'FontName', 'Segoe UI');
        end

            % updating analysis as position changes cla(ax_text);
        axis(ax_text, 'off');

        rectangle('Position', [0.02 0.02 0.96 0.96], 'FaceColor', 'white',
                  ... 'EdgeColor', [0.85 0.85 0.85], 'LineWidth', 1.2, 'Parent',
                  ax_text);

        text(0.50, panel_top, 'STRUCTURAL ANALYSIS', ... 'FontSize', 15,
             'FontWeight', 'bold', 'Parent', ax_text, 'FontName', 'Segoe UI',
             ... 'HorizontalAlignment', 'center', 'Color', [0.12 0.12 0.12]);

        line([panel_left panel_right], [panel_top - 0.03 panel_top - 0.03],
             ... 'Color', [0.75 0.75 0.75], 'LineWidth', 1.2, 'Parent',
             ax_text);

        ypos = panel_top - 0.07;

        % == == == == == == == == == == == == == == == == == == == == == == ==
            == == == == == ==
            = % ADDED : SECTION PROPERTIES BLOCK(constant) % == == == == == ==
              == == == == == == == == == == == == == == == == == == == == == ==
              ==
            = text(panel_left, ypos, 'SECTION PROPERTIES', ... 'FontSize', 12.0,
                   'FontWeight', 'bold', ... 'Parent', ax_text, 'FontName',
                   'Segoe UI', ... 'Color', [0.18 0.18 0.18]);
        ypos = ypos - 0.055;

        text(panel_left + 0.04, ypos, sprintf('\\bar{y} = %.3f mm', ybar),
             ... 'FontSize', 10.5, 'Parent', ax_text, 'FontName', 'Segoe UI',
             ... 'Color', [0.25 0.25 0.25]);
        ypos = ypos - 0.038;

        text(panel_left + 0.04, ypos, sprintf('I = %.3e mm^4', I),
             ... 'FontSize', 10.5, 'Parent', ax_text, 'FontName', 'Segoe UI',
             ... 'Color', [0.25 0.25 0.25]);
        ypos = ypos - 0.050;

        line([panel_left panel_right], [ypos ypos], 'Color', [0.88 0.88 0.88],
             ... 'LineWidth', 0.8, 'Parent', ax_text);
        ypos = ypos - 0.045;
        % == == == == == == == == == == == == == == == == == == == == == == ==
            == == == == == ==
            =

                text(panel_left, ypos, 'CURRENT STATUS', 'FontSize', 12.5,
                     'FontWeight', 'bold', ... 'Parent', ax_text, 'FontName',
                     'Segoe UI', 'Color', [0.18 0.18 0.18]);
        ypos = ypos - 0.06;

        if P_on
          > 0 if FOS_global >= 1.5 fos_color = [0.10 0.58 0.10];
        status_bg = [0.90 0.97 0.90];
        elseif FOS_global >= 1.0 fos_color = [0.88 0.65 0.10];
        status_bg = [1.0 0.97 0.88];
        else fos_color = [0.82 0.08 0.08];
        status_bg = [1.0 0.92 0.92];
        end

            rectangle('Position', [panel_left ypos - 0.012 0.86 0.06],
                      'FaceColor', status_bg, ... 'EdgeColor', 'none', 'Parent',
                      ax_text);

        text(panel_left + 0.04, ypos + 0.018,
             sprintf('Min FOS: %.3f', FOS_global), ... 'FontSize', 14,
             'FontWeight', 'bold', 'Parent', ax_text, ... 'FontName',
             'Segoe UI', 'Color', fos_color);
        ypos = ypos - 0.055;

        text(panel_left + 0.04, ypos,
             sprintf('Critical Mode: %s', critical_mode), ... 'FontSize', 10.2,
             'Parent', ax_text, 'FontName', 'Segoe UI', 'Color',
             [0.25 0.25 0.25]);
        ypos = ypos - 0.042;

        text(panel_left + 0.04, ypos,
             sprintf('Location: x = %.1f mm', fail_location), ... 'FontSize',
             10.2, 'Parent', ax_text, 'FontName', 'Segoe UI', 'Color',
             [0.25 0.25 0.25]);
        else text(panel_left + 0.04, ypos, 'No load currently on bridge',
                  ... 'FontSize', 10.8, 'Color', [0.48 0.48 0.48], 'Parent',
                  ax_text, 'FontName', 'Segoe UI');
        end

            ypos = ypos - 0.06;
        line([panel_left panel_right], [ypos ypos], 'Color', [0.88 0.88 0.88],
             ... 'LineWidth', 0.8, 'Parent', ax_text);
        ypos = ypos - 0.045;

        text(panel_left, ypos, 'LOADING', 'FontSize', 12.5, 'FontWeight',
             'bold', ... 'Parent', ax_text, 'FontName', 'Segoe UI', 'Color',
             [0.18 0.18 0.18]);
        ypos = ypos - 0.055;

        text(panel_left + 0.04, ypos, sprintf('Applied: %.1f N', P_on),
             ... 'FontSize', 10.2, 'Parent', ax_text, 'FontName', 'Segoe UI',
             'Color', [0.25 0.25 0.25]);
        ypos = ypos - 0.04;

        text(panel_left + 0.04, ypos, sprintf('Capacity: %.0f N', total_load),
             ... 'FontSize', 10.2, 'Parent', ax_text, 'FontName', 'Segoe UI',
             'Color', [0.25 0.25 0.25]);

        ypos = ypos - 0.055;
        line([panel_left panel_right], [ypos ypos], 'Color', [0.88 0.88 0.88],
             ... 'LineWidth', 0.8, 'Parent', ax_text);
        ypos = ypos - 0.045;

        text(panel_left, ypos, 'DEFLECTION', 'FontSize', 12.5, 'FontWeight',
             'bold', ... 'Parent', ax_text, 'FontName', 'Segoe UI', 'Color',
             [0.18 0.18 0.18]);
        ypos = ypos - 0.055;

        text(panel_left + 0.04, ypos,
             sprintf('Current: %.3f mm', max_deflection), ... 'FontSize', 10.2,
             'Parent', ax_text, 'FontName', 'Segoe UI', 'Color',
             [0.25 0.25 0.25]);
        ypos = ypos - 0.04;

        text(panel_left + 0.04, ypos,
             sprintf('Span Ratio: L / %.0f', L / max(max_deflection, 0.001)),
             ... 'FontSize', 10.2, 'Parent', ax_text, 'FontName', 'Segoe UI',
             'Color', [0.25 0.25 0.25]);
        ypos = ypos - 0.04;

        text(panel_left + 0.04, ypos,
             sprintf('Historical Max: %.3f mm', max_deflection_history),
             ... 'FontSize', 9.8, 'Parent', ax_text, 'FontName', 'Segoe UI',
             'Color', [0.45 0.45 0.45]);

        ypos = ypos - 0.055;
        line([panel_left panel_right], [ypos ypos], 'Color', [0.88 0.88 0.88],
             ... 'LineWidth', 0.8, 'Parent', ax_text);
        ypos = ypos - 0.045;

        text(panel_left, ypos, 'HISTORICAL MINIMUM', 'FontSize', 12.5,
             'FontWeight', 'bold', ... 'Parent', ax_text, 'FontName',
             'Segoe UI', 'Color', [0.18 0.18 0.18]);
        ypos = ypos - 0.055;

        if min_FOS_history
          < inf text(panel_left + 0.04, ypos,
                     sprintf('Min FOS Ever: %.3f', min_FOS_history),
                     ... 'FontSize', 10.8, 'FontWeight', 'bold', 'Parent',
                     ax_text, ... 'FontName', 'Segoe UI', 'Color',
                     [0.15 0.15 0.68]);
        ypos = ypos - 0.042;

        text(panel_left + 0.04, ypos,
             sprintf('Mode: %s', critical_mode_history), ... 'FontSize', 9.8,
             'Parent', ax_text, 'FontName', 'Segoe UI', 'Color',
             [0.25 0.25 0.25]);
        else text(panel_left + 0.04, ypos, 'Awaiting data...', ... 'FontSize',
                  10.2, 'Color', [0.50 0.50 0.50], 'Parent', ax_text,
                  'FontName', 'Segoe UI');
        end

                %
            == == == == == == == == == == ==
            = % MIN FOS BY MODE section(2 columns) % == == == == == == == == ==
              == ==
            = ypos = ypos - 0.035;
        line([panel_left panel_right], [ypos ypos], 'Color', [0.88 0.88 0.88],
             ... 'LineWidth', 0.8, 'Parent', ax_text);
        ypos = ypos - 0.030;

        text(panel_left, ypos, 'MIN FOS BY MODE (HISTORICAL)', ... 'FontSize',
             11.0, 'FontWeight', 'bold', ... 'Parent', ax_text, 'FontName',
             'Segoe UI', ... 'Color', [0.18 0.18 0.18]);
        ypos = ypos - 0.030;

        % two - column layout nLeft = ceil(nModes / 2);
        % 4 nRight = nModes - nLeft;
        % 3

            x_left = panel_left + 0.04;
        x_right = panel_left + 0.48;
        % shift right column over

                line_step = 0.032;           % vertical spacing per row

    for ii = 1:nModes
        val = min_FOS_mode_history(ii);
        if isfinite (val)
          txt = sprintf('%s: %.3f', modes{ii}, val);
        else
          txt = sprintf('%s: --', modes{ii});
        end

            if val >= 1.5 c = [0.10 0.58 0.10];
        elseif val >= 1.0 c = [0.88 0.65 0.10];
        else c = [0.82 0.08 0.08];
        end

            if ii <= nLeft y_row = ypos - (ii - 1) * line_step;
        x_col = x_left;
        else y_row = ypos - (ii - nLeft - 1) * line_step;
        x_col = x_right;
        end

        text(x_col, y_row, txt, ... 'FontSize', 9.2, 'Parent', ax_text,
             ... 'FontName', 'Segoe UI', 'Color', c);
        end

            drawnow;

        if mod (k, 30)
          == 0 if P_on >
                  0 fprintf(
                      'Frame %d/%d | FOS: %.3f | Load: %.0f N | Deflection: %.2f mm\n',
                      ... k, nSteps, FOS_global, P_on, max_deflection);
        else
          fprintf('Frame %d/%d | No load\n', k, nSteps);
        end end

            % pause once on failure,
            then keep going if failed && ~pause_done pause(2);
        pause_done = true;
        end end

            if ~failed fprintf('\n========== SIMULATION COMPLETE ==========\n');
        fprintf('Bridge survived full train passage!\n');
        fprintf('Minimum FOS achieved: %.3f\n', min_FOS_history);
        fprintf('Critical mode: %s\n', critical_mode_history);
        fprintf('Maximum deflection: %.3f mm (L/%.0f)\n\n',
                ... max_deflection_history,
                L / max(max_deflection_history, 0.001));

        fprintf('--- MIN FOS BY MODE (HISTORICAL) ---\n');
    for
      ii = 1 : nModes if isfinite (min_FOS_mode_history(ii))
                   fprintf('%s: %.3f\n', modes{ii}, min_FOS_mode_history(ii));
    else fprintf('%s: --\n', modes{ii});
    end end fprintf('\n');

    else fprintf('\n========== FINAL FAILURE SUMMARY ==========\n');
    fprintf('Failure location: x = %.1f mm\n', fail_x);
    fprintf('Failure mode: %s\n', fail_mode);
    fprintf('Load at failure: %.1f N (%.1f%% of 400 N)\n', ... P_on_fail,
            100 * P_on_fail / 400);
    fprintf('FOS at failure: %.4f\n', FOS_at_fail);
    fprintf('Maximum deflection reached: %.3f mm (L/%.0f)\n\n',
            ... max_deflection_history, L / max(max_deflection_history, 0.001));

    fprintf('--- MIN FOS BY MODE (HISTORICAL, INCLUDING POST-FAILURE) ---\n');
    for
      ii = 1 : nModes if isfinite (min_FOS_mode_history(ii))
                   fprintf('%s: %.3f\n', modes{ii}, min_FOS_mode_history(ii));
    else fprintf('%s: --\n', modes{ii});
    end end fprintf('\n');
    end

        % % Local funcs function[V, M, RA, RB] =
        compute_SFD_BMD_local(x, L, load_x, load_P) n = numel(x);
    V = zeros(1, n);
    M = zeros(1, n);

    if isempty (load_x)
      RA = 0;
    RB = 0;
    return;
    end

        load_x_snap = interp1(x, x, load_x, 'nearest', 'extrap');
    load_x_snap = max(0, min(L, load_x_snap));

    Ptot = sum(load_P);
    momentA = sum(load_P.*load_x_snap);

    RB = momentA / L;
    RA = Ptot - RB;

    V( :) = RA;
    for
      i = 1 : numel(load_x_snap) V(x >= load_x_snap(i)) =
                  V(x >= load_x_snap(i)) - load_P(i);
    end

        M = cumtrapz(x, V);
    M = M - (M(end) - M(1)) * (x / L) - M(1);
    M(1) = 0;
    M(end) = 0;
    end

        function draw_box_beam_smooth(ax, xline, wline, y_top0, H, zmin, zmax,
                                      faceColor,
                                      edgeColor) yTop = wline + y_top0;
    yBot = yTop - H;
    n = numel(xline);

    [ X_top, Z_top ] = meshgrid(xline, [ zmin, zmax ]);
    Y_top = repmat(yTop, 2, 1);
    surf(ax, X_top, Y_top, Z_top, ... 'FaceColor', faceColor, 'EdgeColor',
         'none', 'FaceAlpha', 0.94);

    [ X_bot, Z_bot ] = meshgrid(xline, [ zmin, zmax ]);
    Y_bot = repmat(yBot, 2, 1);
    surf(ax, X_bot, Y_bot, Z_bot, ... 'FaceColor', faceColor * 0.89,
         'EdgeColor', 'none', 'FaceAlpha', 0.94);

    X_side = repmat(xline, 2, 1);
    Z_left = zmin * ones(2, n);
    Y_left = [yBot; yTop];
    surf(ax, X_side, Y_left, Z_left, ... 'FaceColor', faceColor * 0.82,
         'EdgeColor', 'none', 'FaceAlpha', 0.91);

    Z_right = zmax * ones(2, n);
    Y_right = [yBot; yTop];
    surf(ax, X_side, Y_right, Z_right, ... 'FaceColor', faceColor * 0.82,
         'EdgeColor', 'none', 'FaceAlpha', 0.91);

    xL = xline(1);
    xR = xline(end);

    vertsL = [
      ... xL, yBot(1), zmin; xL, yTop(1), zmin; xL, yTop(1), zmax; xL, yBot(1),
      zmax
    ];
    patch(ax, 'Vertices', vertsL, 'Faces', [1 2 3 4], ... 'FaceColor',
          faceColor * 0.77, 'EdgeColor', edgeColor, ... 'LineWidth', 1.1,
          'FaceAlpha', 0.94);

    vertsR = [
      ... xR, yBot(end), zmin; xR, yTop(end), zmin; xR, yTop(end), zmax;
      xR, yBot(end),
      zmax
    ];
    patch(ax, 'Vertices', vertsR, 'Faces', [1 2 3 4], ... 'FaceColor',
          faceColor * 0.77, 'EdgeColor', edgeColor, ... 'LineWidth', 1.1,
          'FaceAlpha', 0.94);

    plot3(ax, xline, yTop, zmin *ones(1, n), 'Color', edgeColor, 'LineWidth',
          1.4);
    plot3(ax, xline, yTop, zmax *ones(1, n), 'Color', edgeColor, 'LineWidth',
          1.4);
    plot3(ax, xline, yBot, zmin *ones(1, n), 'Color', edgeColor * 0.72,
          'LineWidth', 1.0);
    plot3(ax, xline, yBot, zmax *ones(1, n), 'Color', edgeColor * 0.72,
          'LineWidth', 1.0);

    plot3(ax, [xL xL], [yBot(1) yTop(1)], [zmin zmin], 'Color', edgeColor,
          'LineWidth', 1.1);
    plot3(ax, [xL xL], [yBot(1) yTop(1)], [zmax zmax], 'Color', edgeColor,
          'LineWidth', 1.1);
    plot3(ax, [xR xR], [yBot(end) yTop(end)], [zmin zmin], 'Color', edgeColor,
          'LineWidth', 1.1);
    plot3(ax, [xR xR], [yBot(end) yTop(end)], [zmax zmax], 'Color', edgeColor,
          'LineWidth', 1.1);
    end

        function draw_support_block(ax, x0, height, width) zmin = -width / 2;
    zmax = width / 2;

    verts = [
      x0 - 15, 0, zmin; x0 + 15, 0, zmin; x0 + 15, -height, zmin;
      x0 - 15, -height, zmin; x0 - 15, 0, zmax; x0 + 15, 0, zmax;
      x0 + 15, -height, zmax; x0 - 15, -height,
      zmax
    ];

    faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];

    patch(ax, 'Vertices', verts, 'Faces', faces, ... 'FaceColor',
          [0.38 0.38 0.38], 'EdgeColor', [0.18 0.18 0.18], ... 'LineWidth', 1.3,
          'FaceAlpha', 0.90);
    end
