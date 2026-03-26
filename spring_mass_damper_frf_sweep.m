function spring_mass_damper_frf_sweep()
    % CLEAR WORKSPACE
    clc;  close all;

    %% 1. SYSTEM PARAMETERS
    m = 2.0;       % Mass (kg)
    k = 200.0;     % Spring constant (N/m)
    c = 1.5;       % Damping coefficient (N*s/m)
    
    % Derived Parameters
    wn = sqrt(k/m);          
    zeta = c / (2*sqrt(k*m)); 
    
    % Sweep Parameters
    w_start = 0.5 * wn;      
    w_end   = 2.0 * wn;      
    t_total = 40;            
    ramp_rate = (w_end - w_start) / t_total;
    
    A_base = 0.15; % Base amplitude
    
    % Visual parameters
    L_natural = 1.0;          
    g = 9.81;                 
    delta_static = m * g / k; 
    L_eq = L_natural + delta_static; 
    
    % Initial Conditions
    x0 = [0; 0]; 
    
    %% 2. SOLVE ODE
    ode_fun = @(t, x) get_derivatives(t, x, m, c, k, A_base, w_start, ramp_rate);
    options = odeset('RelTol', 1e-5, 'AbsTol', 1e-5);
    [t_ode, state] = ode45(ode_fun, [0 t_total], x0, options);
    x_mass = state(:, 1);
    
    %% 3. PRE-CALCULATE FRF
    freq_vec = linspace(0, 2.5*wn, 500); 
    r_vec = freq_vec ./ wn;              
    Transmissibility = sqrt( (1 + (2*zeta.*r_vec).^2) ./ ((1 - r_vec.^2).^2 + (2*zeta.*r_vec).^2) );

    %% 4. SETUP FIGURE
    fig = figure('Name', 'Spring-Mass-Damper Sweep', 'Color', 'w', ...
                 'Position', [50, 50, 1200, 700]);
    
    % -- Subplot 1: Animation (Left) --
    ax_anim = subplot(2, 2, [1, 3]);
    hold(ax_anim, 'on'); axis(ax_anim, 'equal');
    ylabel(ax_anim, 'Position (m)');
    xlim(ax_anim, [-1, 1]);
    ylim(ax_anim, [-(L_eq+1.5), 0.5]);
    grid(ax_anim, 'on');
    
    % Static/Dynamic Objects
    base_plate = plot(ax_anim, [-0.6, 0.6], [0, 0], 'k-', 'LineWidth', 5);
    
    % 1. Spring (Left side)
    spring_plot = plot(ax_anim, 0, 0, 'k-', 'LineWidth', 1.5);
    
    % 2. Damper (Right side)
    % We need multiple lines for the damper: Cylinder (3 sides) and Piston (line + head)
    damper_cyl = plot(ax_anim, 0, 0, 'k-', 'LineWidth', 2);
    damper_pis = plot(ax_anim, 0, 0, 'k-', 'LineWidth', 2);
    
    % 3. Mass
    mass_plot = rectangle(ax_anim, 'Position', [0,0,0,0], 'FaceColor', [0.2 0.6 1], 'EdgeColor','k');
    
    txt_freq = text(ax_anim, 0, -L_eq-1.2, '', 'HorizontalAlignment', 'center', 'FontSize', 12);
    
    % -- Subplot 2: Time Response --
    ax_time = subplot(2, 2, 2);
    hold(ax_time, 'on'); grid(ax_time, 'on');
    title(ax_time, 'Time Response');
    xlim(ax_time, [0, t_total]); ylim(ax_time, [-1, 1]);
    trace_mass = plot(ax_time, 0, 0, 'b', 'LineWidth', 1);
    trace_base = plot(ax_time, 0, 0, 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    legend([trace_mass, trace_base], {'Mass', 'Base'});
    
    % -- Subplot 3: FRF --
    ax_frf = subplot(2, 2, 4);
    hold(ax_frf, 'on'); grid(ax_frf, 'on');
    title(ax_frf, 'Transmissibility FRF');
    xlabel(ax_frf, 'Freq (rad/s)'); xlim(ax_frf, [0, max(freq_vec)]);
    plot(ax_frf, freq_vec, Transmissibility, 'k-', 'LineWidth', 1.5);
    xline(ax_frf, wn, '--k', 'Natural Freq');
    frf_dot = plot(ax_frf, 0, 0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    
    %% 5. ANIMATION LOOP
    fps = 30;
    t_anim = 0 : 1/fps : t_total;
    x_anim = interp1(t_ode, x_mass, t_anim);
    
    for i = 1:length(t_anim)
        if ~isvalid(fig), break; end
        
        t_curr = t_anim(i);
        x_curr = x_anim(i);
        
        % Calculate Base Physics
        w_curr = w_start + ramp_rate * t_curr;       
        theta_curr = w_start*t_curr + 0.5*ramp_rate*t_curr^2;
        y_base = A_base * sin(theta_curr);
        
        % Update Animation Elements
        y_mass_abs = -L_eq + x_curr;
        
        % Move Base
        set(base_plate, 'YData', [y_base, y_base]);
        
        % 1. Update Spring (Offset X = -0.2)
        [xs, ys] = get_spring_coords(-0.2, y_base, -0.2, y_mass_abs, 12, 0.15);
        set(spring_plot, 'XData', xs, 'YData', ys);
        
        % 2. Update Damper (Offset X = +0.2)
        % Cylinder attached to Base, Piston attached to Mass
        [xc, yc, xp, yp] = get_damper_coords(0.2, y_base, 0.2, y_mass_abs, 0.15);
        set(damper_cyl, 'XData', xc, 'YData', yc);
        set(damper_pis, 'XData', xp, 'YData', yp);
        
        % 3. Update Mass
        set(mass_plot, 'Position', [-0.4, y_mass_abs-0.2, 0.8, 0.4]); % Wider mass
        set(txt_freq, 'String', sprintf('Freq: %.2f rad/s', w_curr));
        
        % Update Plots
        set(trace_mass, 'XData', t_anim(1:i), 'YData', x_anim(1:i));
        theta_hist = w_start*t_anim(1:i) + 0.5*ramp_rate*t_anim(1:i).^2;
        set(trace_base, 'XData', t_anim(1:i), 'YData', A_base*sin(theta_hist));
        
        r_curr = w_curr / wn;
        Trans_curr = sqrt( (1 + (2*zeta*r_curr)^2) / ((1 - r_curr^2)^2 + (2*zeta*r_curr)^2) );
        set(frf_dot, 'XData', w_curr, 'YData', Trans_curr);
        
        drawnow limitrate;
        pause(1/fps);
    end
end

% --- DERIVATIVES ---
function dxdt = get_derivatives(t, x, m, c, k, A, w_start, ramp)
    w_t = w_start + ramp * t;
    theta_t = w_start*t + 0.5*ramp*t^2;
    y = A * sin(theta_t);
    v_base = A * cos(theta_t) * w_t;
    accel = (c*(v_base - x(2)) + k*(y - x(1))) / m;
    dxdt = [x(2); accel];
end

% --- SPRING COORDS ---
function [x, y] = get_spring_coords(x1, y1, x2, y2, coils, width)
    L = sqrt((x2-x1)^2 + (y2-y1)^2);
    N = coils * 2; 
    s = linspace(0, 1, N+2);
    xs_norm = zeros(size(s));
    xs_norm(2:end-1) = width * (-1).^(1:N);
    theta = atan2(y2-y1, x2-x1) - pi/2;
    xx = xs_norm;
    yy = -linspace(0, -L, length(s)); 
    x = x1 + xx * cos(theta) - yy * sin(theta);
    y = y1 + xx * sin(theta) + yy * cos(theta);
end

% --- DAMPER COORDS ---
function [xc, yc, xp, yp] = get_damper_coords(x1, y1, x2, y2, width)
    % (x1,y1) = Base connection (Cylinder bottom)
    % (x2,y2) = Mass connection (Piston top)
    
    L_tot = sqrt((x2-x1)^2 + (y2-y1)^2);
    theta = atan2(y2-y1, x2-x1) - pi/2;
    
    % Define localized parts
    h_cyl = L_tot * 0.5;  % Cylinder height (halfway)
    h_pis = L_tot * 0.4;  % Piston rod length
    
    % Cylinder: U-shape
    xx_cyl = [-width, -width, width, width];
    yy_cyl = [h_cyl, 0, 0, h_cyl];
    
    % Piston: Line down + horizontal plate
    xx_pis = [0, 0, -width*0.8, width*0.8];
    yy_pis = [0, -h_pis, -h_pis, -h_pis];
    
    % Rotate/Translate Cylinder (Attached to x1,y1)
    xc = x1 + xx_cyl * cos(theta) - yy_cyl * sin(theta);
    yc = y1 + xx_cyl * sin(theta) + yy_cyl * cos(theta);
    
    % Rotate/Translate Piston (Attached to x2,y2)
    % Note: Piston hangs DOWN from x2,y2
    xp = x2 + xx_pis * cos(theta) - yy_pis * sin(theta);
    yp = y2 + xx_pis * sin(theta) + yy_pis * cos(theta);
end