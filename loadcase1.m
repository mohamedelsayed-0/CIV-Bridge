% includes anaylsis similar to the one in python (but not as deep)
% includes a 3D visual of the bridge as it deflects using MAT2
% some boundaries were set to avoid drifting
% for any problems or changes you want to make to the code, you can email me! 
% email: moham.elsayed@mail.utoronto.ca

% note: my way of doing the train load probably wasn't the smartest but it
% works, much better methods can be employed
% this might be quite heavy to run, if you're experiencing issues go to the
% variable nsteps and change it to 180/120 or whatever you need


clear; clc; close all;

% defining section, a new file was made for each different cross-section
% to make it easy to switch from one cross section to another
SEC_FUNC = @section2;
sec = SEC_FUNC();

L       = sec.L;
E       = sec.E;
I       = sec.I;
y_top   = sec.y_top;
y_bot   = sec.y_bot;
Q_web   = sec.Q_web;
Q_glue  = sec.Q_glue;
t_web   = sec.t_web;
t_glue  = sec.t_glue;
caps    = sec.cap;

fprintf('\n========== BRIDGE SIMULATION STARTING ==========\n');
fprintf('Span: %.0f mm | E: %.0f MPa | I: %.1f mm^4\n', L, E, I);
fprintf('Total train load: 400 N (6 axles × %.1f N)\n\n', 400/6);

n  = 600;
x  = linspace(0, L, n+1);

%% Train Load for Loadcase1
total_load = 400;
num_axles  = 6;
P_axle     = total_load / num_axles;

A = sec.spacing_A; % these are the two different spacings between the wheels
B = sec.spacing_B;

offsets = [0, A, A+B, A+B+A, A+B+A+B, A+B+A+B+A];
train_length = offsets(end);

%% setting up visulization of the beam
W_bridge   = 100;
H_bridge   = 60;
zmin = -W_bridge/2;
zmax =  W_bridge/2;

y_top_face0 = 0;

% Animation setup
nSteps   = 260;
front_start = -train_length;
front_end   = L + 60;
front_pos = linspace(front_start, front_end, nSteps);

hFig = figure('Name','Bridge Structural Analysis','NumberTitle','off');
clf(hFig);
set(hFig, 'Color', [0.97 0.97 0.97], 'Position', [50 50 1650 850]);

ax = subplot('Position', [0.04 0.10 0.60 0.85]);
hold(ax,'on'); grid(ax,'on'); box(ax,'on');
view(ax, 130, 25);
camproj(ax,'perspective');
xlabel(ax, 'Longitudinal Position (mm)', 'FontSize', 12, 'FontWeight', 'normal', 'FontName', 'Segoe UI');
ylabel(ax, 'Vertical Deflection (mm)', 'FontSize', 12, 'FontWeight', 'normal', 'FontName', 'Segoe UI');
zlabel(ax, 'Lateral Position (mm)', 'FontSize', 12, 'FontWeight', 'normal', 'FontName', 'Segoe UI');
set(ax, 'FontName', 'Segoe UI', 'FontSize', 10.5, 'Color', [0.99 0.99 0.99]);
set(ax, 'GridAlpha', 0.12, 'GridLineStyle', '-', 'LineWidth', 0.6);

% Side panel for anayalsis
ax_text = subplot('Position', [0.66 0.10 0.32 0.85]);
axis(ax_text,'off');

failed      = false;
fail_x      = NaN;
fail_mode   = '';
FOS_at_fail = NaN;
P_on_fail   = NaN;
eps_small   = 1e-10;

min_FOS_history = inf;
critical_mode_history = '';
max_deflection_history = 0;

for k = 1:nSteps

    % positioning frames
    front_x = front_pos(k);
    axle_x  = front_x + offsets;

    on_span = (axle_x >= 0) & (axle_x <= L);
    load_x  = axle_x(on_span);
    load_P  = P_axle * ones(1, numel(load_x));
    P_on    = sum(load_P);

    % sfd and BMD calcs
    if isempty(load_x)
        V  = zeros(size(x));
        M  = zeros(size(x));
        RA = 0; RB = 0;
    else
        [V, M, RA, RB] = compute_SFD_BMD_local(x, L, load_x, load_P);
    end

    % FOS calcs
    FOS_global    = inf;
    critical_mode = 'None';
    fail_location = NaN;

    if ~failed && P_on > 0
        stress_top = -M .* y_top ./ I;
        stress_bot =  M .* y_bot ./ I;

        stress_bot_tens = max(stress_bot, 0);

        tau_web  = abs(V) .* Q_web  ./ (I * t_web);
        tau_glue = abs(V) .* Q_glue ./ (I * t_glue);

        FOS_tension     = caps.tension          ./ (stress_bot_tens + eps_small);
        FOS_compression = caps.compression      ./ (abs(stress_top) + eps_small);
        FOS_flange_buck = caps.flange_buckling  ./ (abs(stress_top) + eps_small);
        FOS_web_buck    = caps.web_buckling     ./ (abs(stress_top) + eps_small);
        FOS_shear       = caps.shear_board      ./ (tau_web  + eps_small);
        FOS_glue        = caps.shear_glue       ./ (tau_glue + eps_small);
        FOS_shear_buck  = caps.shear_buckling   ./ (tau_web  + eps_small);

        FOS_stack = cat(1, FOS_tension, FOS_compression, FOS_flange_buck, ...
                           FOS_web_buck, FOS_shear, FOS_glue, FOS_shear_buck);

        [FOS_min, mode_idx] = min(FOS_stack, [], 1);
        [FOS_global, iFail] = min(FOS_min);
        fail_location = x(iFail);

        modes = { ...
            'Flexural Tension', ...
            'Flexural Compression', ...
            'Top Flange Buckling', ...
            'Web Buckling', ...
            'Shear in Web', ...
            'Shear in Glue', ...
            'Shear Buckling in Web'};

        critical_mode = modes{mode_idx(iFail)};

        if FOS_global < min_FOS_history
            min_FOS_history = FOS_global;
            critical_mode_history = critical_mode;
        end

        if FOS_global <= 1.0
            failed      = true;
            fail_x      = x(iFail);
            fail_mode   = critical_mode;
            FOS_at_fail = FOS_global;
            P_on_fail   = P_on;

            fprintf('\n!!! BRIDGE FAILURE DETECTED !!!\n');
            fprintf('Frame: %d/%d\n', k, nSteps);
            fprintf('Location: x = %.1f mm\n', fail_x);
            fprintf('Failure mode: %s\n', fail_mode);
            fprintf('Load on span: %.1f N (%.1f%% of total)\n', ...
                    P_on_fail, 100*P_on_fail/total_load);
            fprintf('FOS at failure: %.4f\n\n', FOS_at_fail);
        end
    end

    % Deflection calculation (MAT2)
    w = zeros(size(x));

    if P_on > 0 && max(abs(M)) > 1e-10
        kappa = M ./ (E * I);
        p_fit = polyfit(x, kappa, 1);
        kappa = kappa - polyval(p_fit, x);

        theta = cumtrapz(x, kappa);
        theta = theta - (theta(end) / L) * x;

        w = cumtrapz(x, theta);
        w = w - w(1) - ((w(end) - w(1)) / L) * x;

        w(1) = 0;
        w(end) = 0;

        if max(abs(w)) > L/50
            w = w * (L/50) / max(abs(w));
        end
    end

    max_deflection = max(abs(w));

    if max_deflection > max_deflection_history
        max_deflection_history = max_deflection;
    end

    if max_deflection > 1e-6
        if max_deflection > L/100
            scale_factor = 0.15 * L / max_deflection;
        else
            scale_factor = 0.25 * L / max_deflection;
        end
    else
        scale_factor = 1;
    end

    w_smooth = w * scale_factor;

    % 3D Animation
    cla(ax);
    hold(ax,'on'); grid(ax,'on');

    sup_h = 50;
    draw_support_block(ax, 0, sup_h, W_bridge);
    draw_support_block(ax, L, sup_h, W_bridge);

    if ~failed
        draw_box_beam_smooth(ax, x, w_smooth, y_top_face0, H_bridge, zmin, zmax, ...
                             [0.90 0.95 0.99], [0.28 0.48 0.78]);
    else
        draw_box_beam_smooth(ax, x, w_smooth, y_top_face0, H_bridge, zmin, zmax, ...
                             [1.0 0.72 0.72], [0.87 0.17 0.17]);

        if fail_x >= min(x) && fail_x <= max(x)
            y_fail = interp1(x, w_smooth, fail_x, 'linear', 0) + y_top_face0;
            plot3(ax, fail_x, y_fail, 0, 'rx', 'MarkerSize', 24, 'LineWidth', 5);
            text(fail_x, y_fail + 75, 0, '⚠ FAILURE', ...
                 'FontSize', 13, 'FontWeight', 'bold', 'Color', [0.92 0 0], ...
                 'HorizontalAlignment', 'center', 'Parent', ax, 'FontName', 'Segoe UI');
        end
    end

    if ~isempty(load_x)
        for j = 1:numel(load_x)
            px = load_x(j);

            if px >= min(x) && px <= max(x)
                y_top_here = interp1(x, w_smooth, px, 'linear', 0) + y_top_face0;
            else
                y_top_here = y_top_face0;
            end

            arrow_len = 60;
            quiver3(ax, px, y_top_here + H_bridge + 15, 0, 0, -arrow_len, 0, ...
                    'Color', [0.94 0.14 0.14], 'LineWidth', 3.4, 'MaxHeadSize', 1.6);
        end
    end

    xlim(ax, [-150, L+150]);
    ylim(ax, [-0.4*L, 0.3*L]);
    zlim(ax, [-W_bridge*0.7, W_bridge*0.7]);
    view(ax, 130, 25);

    axis(ax,'vis3d');
    daspect(ax,[1 1 1]);
    camlight(ax,'headlight');
    lighting(ax,'gouraud');

    if ~failed
        title(ax, sprintf('Bridge Structural Model  |  Frame %d/%d', k, nSteps), ...
              'FontWeight', 'normal', 'FontSize', 14, 'FontName', 'Segoe UI', 'Color', [0.15 0.15 0.15]);
    else
        title(ax, sprintf('FAILURE DETECTED  |  x = %.0f mm  |  %s', fail_x, fail_mode), ...
              'FontWeight', 'bold', 'Color', [0.82 0 0], 'FontSize', 14, 'FontName', 'Segoe UI');
    end

    % updating analysis as position changes
    cla(ax_text);
    axis(ax_text,'off');

    rectangle('Position', [0.02 0.02 0.96 0.96], 'FaceColor', 'white', ...
              'EdgeColor', [0.85 0.85 0.85], 'LineWidth', 1.2, 'Parent', ax_text);

    text(0.50, 0.955, 'STRUCTURAL ANALYSIS', ...
         'FontSize', 16, 'FontWeight', 'bold', 'Parent', ax_text, 'FontName', 'Segoe UI', ...
         'HorizontalAlignment', 'center', 'Color', [0.12 0.12 0.12]);
    
    line([0.08 0.92], [0.925 0.925], 'Color', [0.75 0.75 0.75], 'LineWidth', 1.3, 'Parent', ax_text);

        
    ypos = 0.86;
    
    % current Status Section
    text(0.08, ypos, 'CURRENT STATUS', 'FontSize', 13, 'FontWeight', 'bold', ...
         'Parent', ax_text, 'FontName', 'Segoe UI', 'Color', [0.18 0.18 0.18]);
    ypos = ypos - 0.07;

    if P_on > 0
        if FOS_global >= 1.5
            fos_color = [0.10 0.58 0.10];
            status_bg = [0.90 0.97 0.90];
        elseif FOS_global >= 1.0
            fos_color = [0.88 0.65 0.10];
            status_bg = [1.0 0.97 0.88];
        else
            fos_color = [0.82 0.08 0.08];
            status_bg = [1.0 0.92 0.92];
        end
        
        rectangle('Position', [0.08 ypos-0.015 0.84 0.075], 'FaceColor', status_bg, ...
                  'EdgeColor', 'none', 'Parent', ax_text);
        
        text(0.12, ypos+0.025, sprintf('Min FOS: %.3f', FOS_global), ...
             'FontSize', 15, 'FontWeight', 'bold', 'Parent', ax_text, ...
             'FontName', 'Segoe UI', 'Color', fos_color);
        ypos = ypos - 0.065;
        
        text(0.12, ypos, sprintf('Critical Mode: %s', critical_mode), ...
             'FontSize', 10.5, 'Parent', ax_text, 'FontName', 'Segoe UI', 'Color', [0.25 0.25 0.25]);
        ypos = ypos - 0.05;
        
        text(0.12, ypos, sprintf('Location: x = %.1f mm', fail_location), ...
             'FontSize', 10.5, 'Parent', ax_text, 'FontName', 'Segoe UI', 'Color', [0.25 0.25 0.25]);
    else
        text(0.12, ypos, 'No load currently on bridge', ...
             'FontSize', 11.5, 'Color', [0.48 0.48 0.48], 'Parent', ax_text, 'FontName', 'Segoe UI');
    end
    
    ypos = ypos - 0.09;
    line([0.08 0.92], [ypos ypos], 'Color', [0.88 0.88 0.88], 'LineWidth', 0.8, 'Parent', ax_text);
    ypos = ypos - 0.05;
    
    % load info Section
    text(0.08, ypos, 'LOADING', 'FontSize', 13, 'FontWeight', 'bold', ...
         'Parent', ax_text, 'FontName', 'Segoe UI', 'Color', [0.18 0.18 0.18]);
    ypos = ypos - 0.065;
    
    text(0.12, ypos, sprintf('Applied: %.1f N', P_on), ...
         'FontSize', 10.5, 'Parent', ax_text, 'FontName', 'Segoe UI', 'Color', [0.25 0.25 0.25]);
    ypos = ypos - 0.05;
    
    text(0.12, ypos, sprintf('Capacity: %.0f N', total_load), ...
         'FontSize', 10.5, 'Parent', ax_text, 'FontName', 'Segoe UI', 'Color', [0.25 0.25 0.25]);
    
    ypos = ypos - 0.08;
    line([0.08 0.92], [ypos ypos], 'Color', [0.88 0.88 0.88], 'LineWidth', 0.8, 'Parent', ax_text);
    ypos = ypos - 0.05;
    
    % deflection Section
    text(0.08, ypos, 'DEFLECTION', 'FontSize', 13, 'FontWeight', 'bold', ...
         'Parent', ax_text, 'FontName', 'Segoe UI', 'Color', [0.18 0.18 0.18]);
    ypos = ypos - 0.065;
    
    text(0.12, ypos, sprintf('Current: %.3f mm', max_deflection), ...
         'FontSize', 10.5, 'Parent', ax_text, 'FontName', 'Segoe UI', 'Color', [0.25 0.25 0.25]);
    ypos = ypos - 0.05;
    
    text(0.12, ypos, sprintf('Span Ratio: L / %.0f', L/max(max_deflection, 0.001)), ...
         'FontSize', 10.5, 'Parent', ax_text, 'FontName', 'Segoe UI', 'Color', [0.25 0.25 0.25]);
    ypos = ypos - 0.05;
    
    text(0.12, ypos, sprintf('Historical Max: %.3f mm', max_deflection_history), ...
         'FontSize', 10, 'Parent', ax_text, 'FontName', 'Segoe UI', 'Color', [0.45 0.45 0.45]);
    
    ypos = ypos - 0.09;
    line([0.08 0.92], [ypos ypos], 'Color', [0.88 0.88 0.88], 'LineWidth', 0.8, 'Parent', ax_text);
    ypos = ypos - 0.05;
    
    % historical Data section
    text(0.08, ypos, 'HISTORICAL MINIMUM', 'FontSize', 13, 'FontWeight', 'bold', ...
         'Parent', ax_text, 'FontName', 'Segoe UI', 'Color', [0.18 0.18 0.18]);
    ypos = ypos - 0.065;
    
    if min_FOS_history < inf
        hist_color = [0.15 0.15 0.68];
        text(0.12, ypos, sprintf('Min FOS Ever: %.3f', min_FOS_history), ...
             'FontSize', 11.5, 'FontWeight', 'bold', 'Parent', ax_text, ...
             'FontName', 'Segoe UI', 'Color', hist_color);
        ypos = ypos - 0.055;
        
        text(0.12, ypos, sprintf('Mode: %s', critical_mode_history), ...
             'FontSize', 10, 'Parent', ax_text, 'FontName', 'Segoe UI', 'Color', [0.25 0.25 0.25]);
    else
        text(0.12, ypos, 'Awaiting data...', ...
             'FontSize', 11, 'Color', [0.50 0.50 0.50], 'Parent', ax_text, 'FontName', 'Segoe UI');
    end

    if failed
        ypos = 0.11;
        rectangle('Position', [0.08 ypos-0.035 0.84 0.095], 'FaceColor', [1.0 0.92 0.92], ...
                  'EdgeColor', [0.82 0.08 0.08], 'LineWidth', 2.8, 'Parent', ax_text);
        text(0.50, ypos+0.012, '⚠  STRUCTURAL FAILURE DETECTED', ...
             'FontSize', 13.5, 'FontWeight', 'bold', 'Color', [0.82 0 0], ...
             'HorizontalAlignment', 'center', 'Parent', ax_text, 'FontName', 'Segoe UI');
    end

    drawnow;

    if mod(k, 30) == 0 && ~failed
        if P_on > 0
            fprintf('Frame %d/%d | FOS: %.3f | Load: %.0f N | Deflection: %.2f mm\n', ...
                    k, nSteps, FOS_global, P_on, max_deflection);
        else
            fprintf('Frame %d/%d | No load\n', k, nSteps);
        end
    end

    if failed
        pause(2);
        break;
    end
end

if ~failed
    fprintf('\n========== SIMULATION COMPLETE ==========\n');
    fprintf('Bridge survived full train passage!\n');
    fprintf('Minimum FOS achieved: %.3f\n', min_FOS_history);
    fprintf('Critical mode: %s\n', critical_mode_history);
    fprintf('Maximum deflection: %.3f mm (L/%.0f)\n\n', ...
            max_deflection_history, L/max(max_deflection_history, 0.001));
else
    fprintf('\n========== FINAL FAILURE SUMMARY ==========\n');
    fprintf('Failure location: x = %.1f mm\n', fail_x);
    fprintf('Failure mode: %s\n', fail_mode);
    fprintf('Load at failure: %.1f N (%.1f%% of 400 N)\n', ...
            P_on_fail, 100*P_on_fail/400);
    fprintf('FOS at failure: %.4f\n', FOS_at_fail);
    fprintf('Maximum deflection reached: %.3f mm (L/%.0f)\n\n', ...
            max_deflection_history, L/max(max_deflection_history, 0.001));
end

%% Local funcs
function [V, M, RA, RB] = compute_SFD_BMD_local(x, L, load_x, load_P)
    n = numel(x);
    V = zeros(1,n);
    M = zeros(1,n);

    if isempty(load_x)
        RA = 0; RB = 0;
        return;
    end

    load_x_snap = interp1(x, x, load_x, 'nearest', 'extrap');
    load_x_snap = max(0, min(L, load_x_snap));

    Ptot    = sum(load_P);
    momentA = sum(load_P .* load_x_snap);

    RB = momentA / L;
    RA = Ptot - RB;

    V(:) = RA;
    for i = 1:numel(load_x_snap)
        V(x >= load_x_snap(i)) = V(x >= load_x_snap(i)) - load_P(i);
    end

    M = cumtrapz(x, V);
    M = M - (M(end) - M(1)) * (x / L) - M(1);
    M(1)   = 0;
    M(end) = 0;
end

function draw_box_beam_smooth(ax, xline, wline, y_top0, H, zmin, zmax, faceColor, edgeColor)
    yTop = wline + y_top0;
    yBot = yTop - H;
    n    = numel(xline);

    [X_top, Z_top] = meshgrid(xline, [zmin, zmax]);
    Y_top = repmat(yTop, 2, 1);
    surf(ax, X_top, Y_top, Z_top, ...
        'FaceColor', faceColor, 'EdgeColor', 'none', 'FaceAlpha', 0.94);

    [X_bot, Z_bot] = meshgrid(xline, [zmin, zmax]);
    Y_bot = repmat(yBot, 2, 1);
    surf(ax, X_bot, Y_bot, Z_bot, ...
        'FaceColor', faceColor*0.89, 'EdgeColor', 'none', 'FaceAlpha', 0.94);

    X_side = repmat(xline, 2, 1);
    Z_left = zmin * ones(2, n);
    Y_left = [yBot; yTop];
    surf(ax, X_side, Y_left, Z_left, ...
        'FaceColor', faceColor*0.82, 'EdgeColor', 'none', 'FaceAlpha', 0.91);

    Z_right = zmax * ones(2, n);
    Y_right = [yBot; yTop];
    surf(ax, X_side, Y_right, Z_right, ...
        'FaceColor', faceColor*0.82, 'EdgeColor', 'none', 'FaceAlpha', 0.91);

    xL = xline(1);
    xR = xline(end);

    vertsL = [ ...
        xL, yBot(1), zmin;
        xL, yTop(1), zmin;
        xL, yTop(1), zmax;
        xL, yBot(1), zmax ];
    patch(ax, 'Vertices', vertsL, 'Faces', [1 2 3 4], ...
        'FaceColor', faceColor*0.77, 'EdgeColor', edgeColor, ...
        'LineWidth', 1.1, 'FaceAlpha', 0.94);

    vertsR = [ ...
        xR, yBot(end), zmin;
        xR, yTop(end), zmin;
        xR, yTop(end), zmax;
        xR, yBot(end), zmax ];
    patch(ax, 'Vertices', vertsR, 'Faces', [1 2 3 4], ...
        'FaceColor', faceColor*0.77, 'EdgeColor', edgeColor, ...
        'LineWidth', 1.1, 'FaceAlpha', 0.94);

    plot3(ax, xline, yTop, zmin*ones(1,n), 'Color', edgeColor, 'LineWidth', 1.4);
    plot3(ax, xline, yTop, zmax*ones(1,n), 'Color', edgeColor, 'LineWidth', 1.4);
    plot3(ax, xline, yBot, zmin*ones(1,n), 'Color', edgeColor*0.72, 'LineWidth', 1.0);
    plot3(ax, xline, yBot, zmax*ones(1,n), 'Color', edgeColor*0.72, 'LineWidth', 1.0);

    plot3(ax, [xL xL], [yBot(1) yTop(1)], [zmin zmin], 'Color', edgeColor, 'LineWidth', 1.1);
    plot3(ax, [xL xL], [yBot(1) yTop(1)], [zmax zmax], 'Color', edgeColor, 'LineWidth', 1.1);
    plot3(ax, [xR xR], [yBot(end) yTop(end)], [zmin zmin], 'Color', edgeColor, 'LineWidth', 1.1);
    plot3(ax, [xR xR], [yBot(end) yTop(end)], [zmax zmax], 'Color', edgeColor, 'LineWidth', 1.1);
end

function draw_support_block(ax, x0, height, width)
    zmin = -width/2;
    zmax =  width/2;

    verts = [x0-15, 0,       zmin;
             x0+15, 0,       zmin;
             x0+15, -height, zmin;
             x0-15, -height, zmin;
             x0-15, 0,       zmax;
             x0+15, 0,       zmax;
             x0+15, -height, zmax;
             x0-15, -height, zmax];

    faces = [1 2 3 4;
             5 6 7 8;
             1 2 6 5;
             2 3 7 6;
             3 4 8 7;
             4 1 5 8];

    patch(ax, 'Vertices', verts, 'Faces', faces, ...
          'FaceColor', [0.38 0.38 0.38], 'EdgeColor', [0.18 0.18 0.18], ...
          'LineWidth', 1.3, 'FaceAlpha', 0.90);
end
