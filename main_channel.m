%**************************************************************************
% main_channel.m (Presentation Optimized Version)
% Customized main script for Channel Flow analysis
% CFD Lab Assignment 4
%**************************************************************************
close all; clc; clear

% --- 0. Global Figure Settings for Presentation ---
% 统一设置为 Arial 字体，放大字号和线宽，适合 PPT 展示
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 14);
set(0, 'DefaultAxesLineWidth', 1.5);
set(0, 'DefaultLineLineWidth', 2.0);

% --- 1. Initialization and Parameter Setup ---
infilename = 'infile_channel.mat'; 
fprintf('Loading parameters from: %s\n', infilename);

[parm, flow] = build_structs;
[parm, flow] = set_params(parm, flow, infilename);
[parm, flow] = initialize(parm, flow);

% Locate the center of the domain for velocity tracking
i_center = round(parm.m / 2);
j_center = round(parm.n / 2);

% Pre-allocate arrays for convergence tracking
time_history = zeros(parm.ntst, 1);
u_center_history = zeros(parm.ntst, 1);

% Prepare coordinates for plotting
y_coords = linspace(0, parm.yl, parm.n);
[X, Y] = meshgrid(linspace(0, parm.xl, parm.m), linspace(0, parm.yl, parm.n));

% --- Setup for GIF Animation ---
gif_filename = 'channel_flow_development.gif';
fig_anim = figure('Name', 'Flow Development', 'Position', [100, 100, 1000, 450], 'Color', 'w');

fprintf('Starting time integration and recording GIF...\n');

% --- 2. Time Integration Loop ---
for itst = 1 : parm.ntst
    
    % Core Navier-Stokes solver steps
    [flow.rhsu, flow.rhsv] = rhs_ns(parm, flow);
    [flow] = runge_kutta_2d_vec(parm, flow);
    [flow] = direct_press_corr(parm, flow);
    [flow] = project(parm, flow);
    
    current_time = itst * parm.dt;
    time_history(itst) = current_time;
    u_center_history(itst) = flow.u(i_center, j_center);
    
    % --- Capture frames for GIF every 50 steps ---
    if mod(itst, 50) == 0 || itst == 1
        
        % Subplot 1: 2D Velocity Field
        subplot(1, 2, 1);
        contourf(X, Y, flow.u', 30, 'LineStyle', 'none');
        colormap('parula');
        c = colorbar; c.Label.String = 'Velocity U [m/s]'; c.Label.FontSize = 12;
        clim([0, 5.0]); % 固定颜色比例尺防止闪烁
        title(sprintf('2D Velocity Field (t = %.2f s)', current_time), 'FontSize', 16, 'FontWeight', 'bold');
        xlabel('x [m]', 'FontWeight', 'bold'); ylabel('y [m]', 'FontWeight', 'bold');
        axis equal; xlim([0 parm.xl]); ylim([0 parm.yl]);
        
        % Subplot 2: 1D Velocity Profile
        subplot(1, 2, 2);
        plot(flow.u(i_center, :), y_coords, 'b-', 'LineWidth', 2.5);
        title('Velocity Profile Development', 'FontSize', 16, 'FontWeight', 'bold');
        xlabel('Velocity U [m/s]', 'FontWeight', 'bold'); ylabel('Channel Height y [m]', 'FontWeight', 'bold');
        xlim([0, 5.5]); ylim([0, parm.yl]); % 固定坐标轴防止抖动
        grid on;
        
        drawnow;
        
        % Write to GIF
        frame = getframe(fig_anim);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        if itst == 1 || itst == 50
            imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
        else
            imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
        end
    end
    
    if mod(itst, 500) == 0
        fprintf('Progress: Step %d / %d\n', itst, parm.ntst);
    end
end

fprintf('Simulation finished. Generating and saving high-quality static plots...\n');

% --- 3. Post-Processing & Academic Plotting (Static) ---

% Figure 1: Convergence Analysis
fig1 = figure('Name', 'Fig 1: Convergence Analysis', 'Position', [100, 100, 700, 450], 'Color', 'w');
plot(time_history, u_center_history, 'k-', 'LineWidth', 2.5);
title('Convergence Analysis of Centerline Velocity', 'FontSize', 16, 'FontWeight', 'bold');
xlabel('Time [s]', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Centerline Velocity U [m/s]', 'FontSize', 14, 'FontWeight', 'bold');
grid on; set(gca, 'Box', 'on');
hold on; 
xline(3.0, 'r--', 'Steady State Reached', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 2.0, 'FontSize', 12, 'FontWeight', 'bold');
exportgraphics(fig1, 'Fig1_Convergence.png', 'Resolution', 300); % 以300 DPI保存

% Figure 2: Validation of Velocity Profile
u_analytical = (parm.grav / (2 * parm.nu)) .* y_coords .* (parm.yl - y_coords);
u_simulated = flow.u(i_center, :);
fig2 = figure('Name', 'Fig 2: Velocity Profile Validation', 'Position', [150, 150, 600, 650], 'Color', 'w');
plot(u_analytical, y_coords, 'k-', 'LineWidth', 2.5); hold on;
plot(u_simulated, y_coords, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
title('Steady-State Velocity Profile Validation', 'FontSize', 16, 'FontWeight', 'bold');
xlabel('Velocity U [m/s]', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Channel Height y [m]', 'FontSize', 14, 'FontWeight', 'bold');
legend('Analytical Solution', 'Numerical Solution', 'Location', 'northeast', 'FontSize', 12);
grid on; set(gca, 'Box', 'on');
exportgraphics(fig2, 'Fig2_Validation.png', 'Resolution', 300); % 以300 DPI保存

% Figure 3: High-Quality 2D Flow Field with Vectors
fig3 = figure('Name', 'Fig 3: 2D Velocity Field', 'Position', [200, 200, 900, 400], 'Color', 'w');
contourf(X, Y, flow.u', 50, 'LineStyle', 'none'); 
colormap('parula'); 
c = colorbar; c.Label.String = 'Velocity U [m/s]'; c.Label.FontSize = 12;
hold on;
skip = 2; 
quiver(X(1:skip:end, 1:skip:end), Y(1:skip:end, 1:skip:end), ...
       flow.u(1:skip:end, 1:skip:end)', flow.v(1:skip:end, 1:skip:end)', ...
       1.2, 'k', 'LineWidth', 1.2); 
title('Steady-State 2D Velocity Field', 'FontSize', 16, 'FontWeight', 'bold');
xlabel('x [m]', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('y [m]', 'FontSize', 14, 'FontWeight', 'bold');
axis equal; xlim([0 parm.xl]); ylim([0 parm.yl]);
set(gca, 'Box', 'on');
exportgraphics(fig3, 'Fig3_2D_FlowField.png', 'Resolution', 300); % 以300 DPI保存

fprintf('All files successfully saved to your current folder!\n');

% 恢复MATLAB默认字体设置（避免影响你后续的日常使用）
set(0, 'DefaultAxesFontName', 'remove');
set(0, 'DefaultTextFontName', 'remove');
set(0, 'DefaultAxesFontSize', 'remove');
set(0, 'DefaultAxesLineWidth', 'remove');
set(0, 'DefaultLineLineWidth', 'remove');