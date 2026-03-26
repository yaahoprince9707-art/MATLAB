%% MATLAB Code for 2D L-Beam Finite Element Analysis
clear;close all;clc;
% --- 1. Define Parameters and Geometry ---
E = 200e9;       % Young's Modulus (Pa) - Steel
rho = 7850;      % Density (kg/m^3) - Steel
Diameter = 0.05; % Diameter of circular cross-section (m)
L1 = 1.0;        % Length of Column (vertical element, Element 1) (m)
L2 = 0.8;        % Length of Overhang (horizontal element, Element 2) (m)
P = 1000;        % Applied Force at the end of the overhang (N)

% Calculate Cross-sectional Properties for circular section
Area = pi * (Diameter/2)^2; % Area (m^2)
I = pi * (Diameter^4) / 64; % Moment of Inertia (m^4)

% Define Nodes: Total 3 nodes (Node 1=Bottom, Node 2=Corner, Node 3=End)
% Node: [X, Y]
Nodes = [0, 0;  % Node 1 (Bottom, Fixed)
         0, L1; % Node 2 (Intersection/Corner)
         L2, L1]; % Node 3 (Overhang End)

% Define Elements: [Start_Node, End_Node, Length, Angle (rad)]
Elements = [1, 2, L1, pi/2; % Element 1 (Column, vertical, angle 90 deg)
            2, 3, L2, 0];   % Element 2 (Overhang, horizontal, angle 0 deg)

n_nodes = size(Nodes, 1);
n_elem = size(Elements, 1);
DOF_per_node = 3; % u, v, theta
n_dof = n_nodes * DOF_per_node;

% --- 2. Element Stiffness and Mass Matrix Function ---
% Local Beam Stiffness Matrix (Euler-Bernoulli, 6x6)
% DOFs: [u1, v1, theta1, u2, v2, theta2]
k_local = @(L, E, I, A) (E*I/L^3) * ...
    [ (A*L^2/I), 0, 0, -(A*L^2/I), 0, 0;
        0, 12, 6*L, 0, -12, 6*L;
        0, 6*L, 4*L^2, 0, -6*L, 2*L^2;
        -(A*L^2/I), 0, 0, (A*L^2/I), 0, 0;
        0, -12, -6*L, 0, 12, -6*L;
        0, 6*L, 2*L^2, 0, -6*L, 4*L^2 ];

% Local Consistent Mass Matrix (6x6)
m_local = @(L, rho, A) (rho*A*L/420) * ...
    [ 140, 0, 0, 70, 0, 0;
      0, 156, 22*L, 0, 54, -13*L;
      0, 22*L, 4*L^2, 0, 13*L, -3*L^2;
      70, 0, 0, 140, 0, 0;
      0, 54, 13*L, 0, 156, -22*L;
      0, -13*L, -3*L^2, 0, -22*L, 4*L^2 ];

% Rotation Matrix (2D Beam)
T_matrix = @(c, s) ...
    [ c, s, 0, 0, 0, 0;
     -s, c, 0, 0, 0, 0;
      0, 0, 1, 0, 0, 0;
      0, 0, 0, c, s, 0;
      0, 0, 0, -s, c, 0;
      0, 0, 0, 0, 0, 1 ];

% --- 3. Global Assembly ---
K_global = zeros(n_dof);
M_global = zeros(n_dof);

for e = 1:n_elem
    % Element properties
    L = Elements(e, 3);
    angle = Elements(e, 4);
    c = cos(angle);
    s = sin(angle);

    % Calculate local matrices
    k_e = k_local(L, E, I, Area);
    m_e = m_local(L, rho, Area);

    % Rotation Matrix
    T = T_matrix(c, s);

    % Transform local to global
    K_e_global = T' * k_e * T;
    M_e_global = T' * m_e * T;

    % Global DOFs for this element (6 DOFs)
    node_i = Elements(e, 1);
    node_j = Elements(e, 2);
    
    dofs = [ (node_i-1)*DOF_per_node + 1 : node_i*DOF_per_node, ...
             (node_j-1)*DOF_per_node + 1 : node_j*DOF_per_node ];
    
    % Assembly
    K_global(dofs, dofs) = K_global(dofs, dofs) + K_e_global;
    M_global(dofs, dofs) = M_global(dofs, dofs) + M_e_global;
end

% Global Force Vector F
F_global = zeros(n_dof, 1);
% Force P acts vertically (-Y direction) at Node 3 (DOFs: u3, v3, theta3)
F_global(n_dof - 1) = -P; % v3 is the last-1 DOF

% --- 4. Apply Boundary Conditions and Solve ---
% Fixed Support at Node 1 (Bottom)
% DOFs to restrain: u1, v1, theta1 (DOFs 1, 2, 3)
prescribed_dofs = [1, 2, 3]; 

% Free DOFs
all_dofs = 1:n_dof;
free_dofs = setdiff(all_dofs, prescribed_dofs);

% Reduced System: K_ff * D_f = F_f
K_ff = K_global(free_dofs, free_dofs);
F_f = F_global(free_dofs);

% Solve for the unknown displacements
D_f = K_ff \ F_f;

% Assemble the full displacement vector
D_global = zeros(n_dof, 1);
D_global(free_dofs) = D_f;

% Extract nodal displacements [u1, v1, theta1, u2, v2, theta2, u3, v3, theta3]
D_nodes = reshape(D_global, DOF_per_node, n_nodes)';
D_u = D_nodes(:, 1); % u displacement
D_v = D_nodes(:, 2); % v displacement
D_theta = D_nodes(:, 3); % theta rotation

% --- 5. Calculate Internal Forces in Each Element ---
fprintf('\n--- Internal Element Forces ---\n');
% The reaction forces/moments at the local nodes (P1, V1, M1, P2, V2, M2)
Internal_Forces = zeros(n_elem, 6);

for e = 1:n_elem
    L = Elements(e, 3);
    angle = Elements(e, 4);
    c = cos(angle);
    s = sin(angle);

    % Global DOFs for this element (6 DOFs)
    node_i = Elements(e, 1);
    node_j = Elements(e, 2);
    dofs = [ (node_i-1)*DOF_per_node + 1 : node_i*DOF_per_node, ...
             (node_j-1)*DOF_per_node + 1 : node_j*DOF_per_node ];

    % Element displacement vector (global)
    d_global = D_global(dofs);

    % Transform global displacement to local
    T = T_matrix(c, s);
    d_local = T * d_global;

    % Calculate local internal forces: f_local = k_local * d_local
    k_e = k_local(L, E, I, Area);
    f_local = k_e * d_local;

    Internal_Forces(e, :) = f_local;

    % Print results (P=Axial Force, V=Shear Force, M=Bending Moment)
    fprintf('Element %d (Nodes %d -> %d):\n', e, node_i, node_j);
    fprintf('  Node 1: Axial=%.2f N, Shear=%.2f N, Moment=%.2f Nm\n', f_local(1), f_local(2), f_local(3));
    fprintf('  Node 2: Axial=%.2f N, Shear=%.2f N, Moment=%.2f Nm\n', f_local(4), f_local(5), f_local(6));
end

% --- 6. Plot Deformed and Undeformed Structure ---
Deformation_Scale = 100; % Scale up displacement for visualization
New_X = Nodes(:, 1) + D_u * Deformation_Scale;
New_Y = Nodes(:, 2) + D_v * Deformation_Scale;

figure;
hold on;
title('2D L-Beam Deformed vs. Undeformed Shape');
xlabel('X Position (m)');
ylabel('Y Position (m)');

% Plot Undeformed (black dashed)
plot(Nodes(:, 1), Nodes(:, 2), 'k--', 'LineWidth', 2, 'DisplayName', 'Undeformed');
scatter(Nodes(:, 1), Nodes(:, 2), 50, 'k', 'filled'); % Plot nodes

% Plot Deformed (red solid)
plot(New_X, New_Y, 'r-', 'LineWidth', 2, 'DisplayName', 'Deformed (Scaled)');
scatter(New_X, New_Y, 50, 'r', 'filled'); % Plot deformed nodes

% Annotate Force P at Node 3
text(Nodes(3, 1), Nodes(3, 2), sprintf('  Force P = %.0f N', P), 'VerticalAlignment', 'bottom');

% Annotate Fixed Support at Node 1
text(Nodes(1, 1), Nodes(1, 2), 'Fixed Support', 'HorizontalAlignment', 'right');

grid on;
axis equal;
legend('Location', 'best');
hold off;