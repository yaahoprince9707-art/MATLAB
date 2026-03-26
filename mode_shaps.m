%% beam_modes_udl.m
% FEM modal analysis (Euler-Bernoulli beam) + static deflection under UDL
clear; close all; clc;

%% Parameters
L = 2.0;               % total length (m)
ne = 40;               % number of elements (increase for higher accuracy)
nn = ne + 1;           % number of nodes
le = L / ne;           % element length

E = 210e9;             % Young's modulus (Pa)
rho = 7850;            % density (kg/m^3)
b = 0.02;              % width (m)
h = 0.02;              % thickness (m)
A = b * h;             % cross-sectional area (m^2)
I = (b * h^3) / 12;    % second moment of area (m^4)

q = 1000;              % UDL (N/m) acting downward (positive number) 

%% DOF setup
ndof = 2 * nn;         % DOFs: [w1 theta1 w2 theta2 ...]
K = zeros(ndof);
M = zeros(ndof);
F = zeros(ndof,1);     % global nodal force vector (for UDL)

%% Element matrices (consistent)
for e = 1:ne
    % element DOF indices
    n1 = e;
    n2 = e+1;
    dofs = [2*n1-1, 2*n1, 2*n2-1, 2*n2]; % [w1 th1 w2 th2]
    
    % stiffness matrix (4x4) for Euler-Bernoulli beam element
    ke = E*I / le^3 * [  12,   6*le,  -12,   6*le;
                         6*le, 4*le^2, -6*le, 2*le^2;
                        -12,  -6*le,   12,  -6*le;
                         6*le, 2*le^2, -6*le, 4*le^2 ];
                     
    % consistent mass matrix (4x4)
    me = rho*A*le/420 * [156, 22*le, 54, -13*le;
                         22*le, 4*le^2, 13*le, -3*le^2;
                         54, 13*le, 156, -22*le;
                        -13*le, -3*le^2, -22*le, 4*le^2 ];
    
    % consistent nodal load vector for UDL q per length (downwards)
    % standard element consistent load: fe = q*le/12 * [6; le; 6; -le]
    fe = q * le / 12 * [6; le; 6; -le];
    % note: sign convention - we'll treat downward positive and plot accordingly
    
    % assemble
    K(dofs,dofs) = K(dofs,dofs) + ke;
    M(dofs,dofs) = M(dofs,dofs) + me;
    F(dofs) = F(dofs) + fe;
end

%% Apply boundary conditions (simply supported)
% For simple support at node 1 and node nn:
% - transverse displacement w = 0 at node 1 and node nn  -> fix DOFs 1 and 2*nn-1
% - rotations are free (do NOT fix rotation DOFs)
fixedDOF = [1, 2*nn-1];

allDOF = 1:ndof;
freeDOF = setdiff(allDOF, fixedDOF);

%% Static solution under UDL (K * u = F)
% We remove fixed displacement DOFs and solve for the rest.
Kff = K(freeDOF, freeDOF);
Ff = F(freeDOF);

uf = Kff \ Ff;   % reduced solution
u = zeros(ndof,1);
u(freeDOF) = uf;

% Extract transverse displacements only for plotting static deflection
w_static = u(1:2:end);

%% Modal analysis: solve generalized eigenvalue problem K * phi = omega^2 * M * phi
Krr = K(freeDOF, freeDOF);
Mrr = M(freeDOF, freeDOF);

% Solve eigenproblem (compute first few modes)
nModes = 6; % how many modes to compute / plot
opts.maxit = 500;
% Use eigs for efficiency for larger matrices
[phi_r, D] = eigs(Krr, Mrr, nModes, 'SM'); % 'SM' -> smallest magnitude eigenvalues (lowest freq)
omega2 = diag(D);
% remove any small negative numerical noise
omega2 = real(omega2);
omega = sqrt(omega2);
freq = omega / (2*pi);

% build full-mode vectors (including fixed DOFs = 0 for w at supports)
Phi = zeros(ndof, nModes);
for i = 1:nModes
    Phi(freeDOF,i) = real(phi_r(:,i));
end

% extract transverse DOFs for plotting
xnodes = linspace(0, L, nn)';
modes_w = Phi(1:2:end, :);

% normalize modes: max amplitude = 1 (or sign such that first entry positive)
for i = 1:nModes
    modes_w(:,i) = modes_w(:,i) / max(abs(modes_w(:,i)));
    % enforce a consistent sign (make average positive)
    if mean(modes_w(:,i)) < 0
        modes_w(:,i) = -modes_w(:,i);
    end
end

%% Plot static deflection and first nModes mode shapes
figure('Position',[100 100 1000 700]);
nm_to_show = min(nModes, 6);
for i = 1:nm_to_show
    subplot(3,2,i);
    plot(xnodes, modes_w(:,i), '-o','LineWidth',1.2,'MarkerSize',4);
    hold on;
    % overlay static deflection scaled for visibility
    scaleStatic = 5 * max(abs(modes_w(:,i))) / max(abs(w_static)+eps); % scale factor for visualization
    plot(xnodes, -scaleStatic * w_static, '--','LineWidth',1.1); % negative sign to show downward deflection consistent
    hold off;
    title(sprintf('Mode %d — f = %.3f Hz', i, freq(i)));
    xlabel('x (m)');
    ylabel('Normalized transverse displacement');
    legend('Mode shape (FEM eigenvector)','Scaled static deflection (UDL)','Location','best');
    grid on;
end
sgtitle('Simply Supported Beam — FEM Mode Shapes and Scaled Static Deflection (UDL)');

%% Optionally: print frequencies
fprintf('Natural frequencies (Hz) for first %d modes:\n', nm_to_show);
for i = 1:nm_to_show
    fprintf(' Mode %d: %.6f Hz\n', i, freq(i));
end
