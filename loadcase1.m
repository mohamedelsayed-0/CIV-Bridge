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
x  = linspace(0, L, 600);

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
set(hFig, 'Color', 'w', 'Position', [100 100 1600 800]);

% Main 3D view
ax = subplot('Position', [0.05 0.12 0.58 0.82]);
hold(ax,'on'); grid(ax,'on'); box(ax,'on');
view(ax, 130, 25);
camproj(ax,'perspective');
xlabel(ax, 'Longitudinal Position (mm)', 'FontSize', 11, 'FontWeight', 'normal', 'FontName', 'Arial');
ylabel(ax, 'Vertical Deflection (mm)', 'FontSize', 11, 'FontWeight', 'normal', 'FontName', 'Arial');
zlabel(ax, 'Lateral Position (mm)', 'FontSize', 11, 'FontWeight', 'normal', 'FontName', 'Arial');
set(ax, 'FontName', 'Arial', 'FontSize', 10);

% Side panel for anayalsis
ax_text = subplot('Position', [0.67 0.12 0.30 0.82]);
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
        S_top = abs(M) .* y_top ./ I;
        S_bot = abs(M) .* y_bot ./ I;

        tau_web  = abs(V) .* Q_web  ./ (I * t_web);
        tau_glue = abs(V) .* Q_glue ./ (I * t_glue);

        FOS_tension     = caps.tension          ./ (S_bot    + eps_small);
        FOS_compression = caps.compression      ./ (S_top    + eps_small);
        FOS_flange_buck = caps.flange_buckling  ./ (S_top    + eps_small);
        FOS_web_buck    = caps.web_buckling     ./ (S_top    + eps_small);
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
    set(ax, 'GridAlpha', 0.15, 'GridLineStyle', '-', 'LineWidth', 0.5);

    sup_h = 50;
    draw_support_block(ax, 0, sup_h, W_bridge);
    draw_support_block(ax, L, sup_h, W_bridge);

    if ~failed
        draw_box_beam_smooth(ax, x, w_smooth, y_top_face0, H_bridge, zmin, zmax, ...
                             [0.88 0.94 0.99], [0.25 0.45 0.75]);
    else
        draw_box_beam_smooth(ax, x, w_smooth, y_top_face0, H_bridge, zmin, zmax, ...
                             [1.0 0.70 0.70], [0.85 0.15 0.15]);

        if fail_x >= min(x) && fail_x <= max(x)
            y_fail = interp1(x, w_smooth, fail_x, 'linear', 0) + y_top_face0;
            plot3(ax, fail_x, y_fail, 0, 'rx', 'MarkerSize', 22, 'LineWidth', 4.5);
            text(fail_x, y_fail + 70, 0, '⚠ FAILURE', ...
                 'FontSize', 12, 'FontWeight', 'bold', 'Color', [0.9 0 0], ...
                 'HorizontalAlignment', 'center', 'Parent', ax, 'FontName', 'Arial');
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
                    'Color', [0.92 0.12 0.12], 'LineWidth', 3.2, 'MaxHeadSize', 1.5);
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
        title(ax, sprintf('Bridge Structural Model | Frame %d/%d', k, nSteps), ...
              'FontWeight', 'bold', 'FontSize', 13, 'FontName', 'Arial');
    else
        title(ax, sprintf('FAILURE DETECTED | x = %.0f mm | %s', fail_x, fail_mode), ...
              'FontWeight', 'bold', 'Color', [0.8 0 0], 'FontSize', 13, 'FontName', 'Arial');
    end

    % updating analysis as position changes
    cla(ax_text);
    axis(ax_text,'off');

    % title
    text(0.05, 0.96, 'STRUCTURAL ANALYSIS', ...
         'FontSize', 15, 'FontWeight', 'bold', 'Parent', ax_text, 'FontName', 'Arial');
    
    % divider line
    line([0.05 0.95], [0.935 0.935], 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5, 'Parent', ax_text);

    ypos = 0.88;
    
    % current Status Section
    text(0.05, ypos, 'CURRENT STATUS', 'FontSize', 12, 'FontWeight', 'bold', ...
         'Parent', ax_text, 'FontName', 'Arial', 'Color', [0.2 0.2 0.2]);
    ypos = ypos - 0.06;

    if P_on > 0
        if FOS_global >= 1.5
            fos_color = [0.13 0.55 0.13];  % green
        elseif FOS_global >= 1.0
            fos_color = [0.85 0.65 0.13];  % amber
        else
            fos_color = [0.8 0.1 0.1];     % red
        end
        
        text(0.08, ypos, sprintf('Min FOS: %.3f', FOS_global), ...
             'FontSize', 13, 'FontWeight', 'bold', 'Parent', ax_text, ...
             'FontName', 'Arial', 'Color', fos_color);
        ypos = ypos - 0.05;
        
        text(0.08, ypos, sprintf('Mode: %s', critical_mode), ...
             'FontSize', 10, 'Parent', ax_text, 'FontName', 'Arial', 'Color', [0.3 0.3 0.3]);
        ypos = ypos - 0.045;
        
        text(0.08, ypos, sprintf('Location: x = %.1f mm', fail_location), ...
             'FontSize', 10, 'Parent', ax_text, 'FontName', 'Arial', 'Color', [0.3 0.3 0.3]);
    else
        text(0.08, ypos, 'No load on bridge', ...
             'FontSize', 11, 'Color', [0.5 0.5 0.5], 'Parent', ax_text, 'FontName', 'Arial');
    end
    
    ypos = ypos - 0.08;
    
    % load info Section
    text(0.05, ypos, 'LOADING', 'FontSize', 12, 'FontWeight', 'bold', ...
         'Parent', ax_text, 'FontName', 'Arial', 'Color', [0.2 0.2 0.2]);
    ypos = ypos - 0.06;
    
    text(0.08, ypos, sprintf('Applied: %.1f N / %.0f N', P_on, total_load), ...
         'FontSize', 10, 'Parent', ax_text, 'FontName', 'Arial');
    ypos = ypos - 0.045;
    
    text(0.08, ypos, sprintf('Progress: %.1f%%', 100*P_on/total_load), ...
         'FontSize', 10, 'Parent', ax_text, 'FontName', 'Arial');
    
    ypos = ypos - 0.08;
    
    % deflection Section
    text(0.05, ypos, 'DEFLECTION', 'FontSize', 12, 'FontWeight', 'bold', ...
         'Parent', ax_text, 'FontName', 'Arial', 'Color', [0.2 0.2 0.2]);
    ypos = ypos - 0.06;
    
    text(0.08, ypos, sprintf('Max: %.3f mm', max_deflection), ...
         'FontSize', 10, 'Parent', ax_text, 'FontName', 'Arial');
    ypos = ypos - 0.045;
    
    text(0.08, ypos, sprintf('Span Ratio: L/%.0f', L/max(max_deflection, 0.001)), ...
         'FontSize', 10, 'Parent', ax_text, 'FontName', 'Arial');
    ypos = ypos - 0.045;
    
    text(0.08, ypos, sprintf('Historical Max: %.3f mm', max_deflection_history), ...
         'FontSize', 9, 'Parent', ax_text, 'FontName', 'Arial', 'Color', [0.4 0.4 0.4]);
    
    ypos = ypos - 0.10;
    
    % historical Data section
    line([0.05 0.95], [ypos+0.02 ypos+0.02], 'Color', [0.7 0.7 0.7], 'LineWidth', 1, 'Parent', ax_text);
    
    text(0.05, ypos, 'HISTORICAL MINIMUM', 'FontSize', 12, 'FontWeight', 'bold', ...
         'Parent', ax_text, 'FontName', 'Arial', 'Color', [0.2 0.2 0.2]);
    ypos = ypos - 0.06;
    
    if min_FOS_history < inf
        hist_color = [0.2 0.2 0.7];
        text(0.08, ypos, sprintf('FOS: %.3f', min_FOS_history), ...
             'FontSize', 12, 'FontWeight', 'bold', 'Parent', ax_text, ...
             'FontName', 'Arial', 'Color', hist_color);
        ypos = ypos - 0.05;
        
        text(0.08, ypos, sprintf('Mode: %s', critical_mode_history), ...
             'FontSize', 9, 'Parent', ax_text, 'FontName', 'Arial', 'Color', [0.3 0.3 0.3]);
    else
        text(0.08, ypos, 'Awaiting data...', ...
             'FontSize', 10, 'Color', [0.5 0.5 0.5], 'Parent', ax_text, 'FontName', 'Arial');
    end

    % failure msg
    if failed
        ypos = 0.08;
        rectangle('Position', [0.05 ypos-0.04 0.90 0.10], 'FaceColor', [1 0.9 0.9], ...
                  'EdgeColor', [0.8 0.1 0.1], 'LineWidth', 2.5, 'Parent', ax_text);
        text(0.50, ypos+0.01, '⚠ STRUCTURAL FAILURE', ...
             'FontSize', 13, 'FontWeight', 'bold', 'Color', [0.8 0 0], ...
             'HorizontalAlignment', 'center', 'Parent', ax_text, 'FontName', 'Arial');
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
        'FaceColor', faceColor, 'EdgeColor', 'none', 'FaceAlpha', 0.92);

    [X_bot, Z_bot] = meshgrid(xline, [zmin, zmax]);
    Y_bot = repmat(yBot, 2, 1);
    surf(ax, X_bot, Y_bot, Z_bot, ...
        'FaceColor', faceColor*0.88, 'EdgeColor', 'none', 'FaceAlpha', 0.92);

    X_side = repmat(xline, 2, 1);
    Z_left = zmin * ones(2, n);
    Y_left = [yBot; yTop];
    surf(ax, X_side, Y_left, Z_left, ...
        'FaceColor', faceColor*0.80, 'EdgeColor', 'none', 'FaceAlpha', 0.88);

    Z_right = zmax * ones(2, n);
    Y_right = [yBot; yTop];
    surf(ax, X_side, Y_right, Z_right, ...
        'FaceColor', faceColor*0.80, 'EdgeColor', 'none', 'FaceAlpha', 0.88);

    xL = xline(1);
    xR = xline(end);

    vertsL = [ ...
        xL, yBot(1), zmin;
        xL, yTop(1), zmin;
        xL, yTop(1), zmax;
        xL, yBot(1), zmax ];
    patch(ax, 'Vertices', vertsL, 'Faces', [1 2 3 4], ...
        'FaceColor', faceColor*0.75, 'EdgeColor', edgeColor, ...
        'LineWidth', 1.0, 'FaceAlpha', 0.92);

    vertsR = [ ...
        xR, yBot(end), zmin;
        xR, yTop(end), zmin;
        xR, yTop(end), zmax;
        xR, yBot(end), zmax ];
    patch(ax, 'Vertices', vertsR, 'Faces', [1 2 3 4], ...
        'FaceColor', faceColor*0.75, 'EdgeColor', edgeColor, ...
        'LineWidth', 1.0, 'FaceAlpha', 0.92);

    plot3(ax, xline, yTop, zmin*ones(1,n), 'Color', edgeColor, 'LineWidth', 1.3);
    plot3(ax, xline, yTop, zmax*ones(1,n), 'Color', edgeColor, 'LineWidth', 1.3);
    plot3(ax, xline, yBot, zmin*ones(1,n), 'Color', edgeColor*0.7, 'LineWidth', 0.9);
    plot3(ax, xline, yBot, zmax*ones(1,n), 'Color', edgeColor*0.7, 'LineWidth', 0.9);

    plot3(ax, [xL xL], [yBot(1) yTop(1)], [zmin zmin], 'Color', edgeColor, 'LineWidth', 1.0);
    plot3(ax, [xL xL], [yBot(1) yTop(1)], [zmax zmax], 'Color', edgeColor, 'LineWidth', 1.0);
    plot3(ax, [xR xR], [yBot(end) yTop(end)], [zmin zmin], 'Color', edgeColor, 'LineWidth', 1.0);
    plot3(ax, [xR xR], [yBot(end) yTop(end)], [zmax zmax], 'Color', edgeColor, 'LineWidth', 1.0);
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
          'FaceColor', [0.35 0.35 0.35], 'EdgeColor', [0.15 0.15 0.15], ...
          'LineWidth', 1.2, 'FaceAlpha', 0.88);
end