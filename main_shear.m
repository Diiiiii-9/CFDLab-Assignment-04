%**************************************************************************
% main_shear.m
% Customized main script for Unsteady Shear Layer analysis
% CFD Lab Assignment 4
%**************************************************************************
close all; clc; clear

% --- 0. Global Plotting Settings (Academic Quality) ---
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 14);
set(0, 'DefaultAxesLineWidth', 1.5);
set(0, 'DefaultLineLineWidth', 2.0);

% --- 1. Initialization and Parameter Setup ---
infilename = 'infile_shear.mat'; 
fprintf('Loading parameters from: %s\n', infilename);

[parm, flow] = build_structs;
[parm, flow] = set_params(parm, flow, infilename);
[parm, flow] = initialize(parm, flow);

% Pre-allocate array for Global Kinetic Energy tracking
time_history = zeros(parm.ntst, 1);
KE_history = zeros(parm.ntst, 1);

% Get coordinates for plotting 
[X, Y] = meshgrid(linspace(0, parm.xl, parm.m), linspace(0, parm.yl, parm.n));
dx = parm.xl / parm.m;
dy = parm.yl / parm.n;

% --- Setup for GIF Animation ---
fprintf('Setting up GIF generation for Vorticity Field...\n');
fig_gif = figure('Name', 'Vorticity Roll-up Animation', 'Position', [100, 200, 700, 600], 'Color', 'w');
gif_filename = 'shear_layer_vorticity.gif';
if exist(gif_filename, 'file')
    delete(gif_filename); 
end

fprintf('Starting time integration...\n');

% --- 2. Time Integration Loop ---
for itst = 1 : parm.ntst
    
    % Core Navier-Stokes solver steps
    [flow.rhsu, flow.rhsv] = rhs_ns(parm, flow);
    [flow] = runge_kutta_2d_vec(parm, flow);
    [flow] = direct_press_corr(parm, flow);
    [flow] = project(parm, flow);
    
    % Record time and Total Kinetic Energy (E_k = 0.5 * sum(u^2 + v^2))
    time_history(itst) = itst * parm.dt;
    KE_history(itst) = 0.5 * sum(sum(flow.u.^2 + flow.v.^2)) * dx * dy;
    
    % --- GIF Frame Capture: Vorticity Evolution (every 50 steps) ---
    if mod(itst, 50) == 0
        figure(fig_gif); clf;
        
        % Calculate Vorticity using MATLAB's built-in curl function
        [curlz, ~] = curl(X, Y, flow.u', flow.v');
        
        % Plot Vorticity Field
        contourf(X, Y, curlz, 40, 'LineStyle', 'none');
        colormap(parula); % Use a diverging colormap (red/blue) for positive/negative vorticity if available, else parula
        c = colorbar; c.Label.String = 'Vorticity \omega [1/s]';
        caxis([-5 5]); % Fix color scale to prevent flashing (adjust limits based on your data)
        
        title(sprintf('Vorticity Field (Time = %.2f s)', itst*parm.dt), 'FontSize', 16, 'FontWeight', 'bold');
        xlabel('x [m]'); ylabel('y [m]'); axis equal; 
        xlim([0 parm.xl]); ylim([0 parm.yl]);
        set(gca, 'LineWidth', 1.2, 'Box', 'on');
        drawnow;
        
        % Append to GIF
        frame = getframe(fig_gif);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        if itst == 50
            imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
        else
            imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
        end
    end
    
    if mod(itst, 500) == 0
        fprintf('Progress: Step %d / %d\n', itst, parm.ntst);
    end
end

close(fig_gif); 
fprintf('Simulation finished. Generating high-quality academic plots...\n');

% --- 3. Post-Processing & Academic Plotting ---

% -------------------------------------------------------------------------
% Figure 1: Total Kinetic Energy Decay 
% -------------------------------------------------------------------------
fig1 = figure('Name', 'Fig 1: Kinetic Energy Decay', 'Position', [100, 100, 600, 450], 'Color', 'w');
plot(time_history, KE_history / KE_history(1), 'k-', 'LineWidth', 2.0); % Normalized
title('Decay of Normalized Total Kinetic Energy', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Time [s]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('E_k(t) / E_k(0)', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'Box', 'on');
exportgraphics(fig1, 'Fig1_KineticEnergy_Shear.png', 'Resolution', 300);

% -------------------------------------------------------------------------
% Figure 2: Final Vorticity Field (The "Gold Standard" Plot)
% -------------------------------------------------------------------------
[curlz, ~] = curl(X, Y, flow.u', flow.v');

fig2 = figure('Name', 'Fig 2: Final Vorticity Field', 'Position', [750, 100, 600, 500], 'Color', 'w');
contourf(X, Y, curlz, 50, 'LineStyle', 'none'); 
colormap('parula'); 
c = colorbar; c.Label.String = 'Vorticity \omega [1/s]';
title('Final Vorticity Field (Vortex Roll-up)', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('x [m]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('y [m]', 'FontSize', 12, 'FontWeight', 'bold');
axis equal; xlim([0 parm.xl]); ylim([0 parm.yl]);
set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'Box', 'on');
exportgraphics(fig2, 'Fig2_VorticityField_Shear.png', 'Resolution', 300);

% -------------------------------------------------------------------------
% Figure 3: Final Streamlines and Velocity Magnitude
% -------------------------------------------------------------------------
velocity_mag = sqrt(flow.u'.^2 + flow.v'.^2);

fig3 = figure('Name', 'Fig 3: Final Streamlines', 'Position', [200, 550, 600, 500], 'Color', 'w');
% Plot velocity magnitude as background
imagesc([0 parm.xl], [0 parm.yl], velocity_mag); 
set(gca, 'YDir', 'normal'); % Correct axis direction for imagesc
colormap(gca, 'parula');
c = colorbar; c.Label.String = 'Velocity Magnitude |V| [m/s]';
hold on;

% Overlay streamlines to show vortex topology
streamslice(X, Y, flow.u', flow.v', 2); % '2' controls streamline density
set(findobj(gca, 'Type', 'line'), 'Color', 'w', 'LineWidth', 1.0); % White streamlines

title('Final Streamlines and Velocity Topology', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('x [m]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('y [m]', 'FontSize', 12, 'FontWeight', 'bold');
axis equal; xlim([0 parm.xl]); ylim([0 parm.yl]);
set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'Box', 'on');
exportgraphics(fig3, 'Fig3_Streamlines_Shear.png', 'Resolution', 300);

fprintf('Success: High-quality PNG images and Vorticity GIF have been saved.\n');