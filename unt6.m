
clear; close all; clc; format short g

m = [1; 1; 1]; 
k1 = 1000; k2 = 1000; k3 = 1000; k4 = 1000;

mc = 1; kc = 1000;%chain

%Beam geometry
L = 1.0;            % Length (m)
E = 210e9;          % Young's modulus (Pa)
rho = 7800;         % density (kg/m^3)
b = 0.02; h = 0.01; % cross-section (m)
A = b*h;
I = b*h^3/12;

% mode count to extract
nModes = 6; % compute more modes then pick first 3 for convergence

%% PART A: Discrete 3-DOF system
K3 = [k1+k2,   -k2,     0;
      -k2,   k2+k3,   -k3;
       0,     -k3,   k3+k4];
M3 = diag(m);

[phi3, D3] = eig(K3, M3);
omega2 = diag(D3);
[omega2_sorted, idx] = sort(omega2);
phi3 = phi3(:, idx);
omega3 = sqrt(omega2_sorted);
freq3 = omega3/(2*pi);

% fprintf('3-DOF discrete natural frequencies (Hz):\n');
% disp(freq3(1:3))

% Normalize mode shapes (by max amplitude for plotting)
for i=1:3
    phi3(:,i) = phi3(:,i)/max(abs(phi3(:,i)));
end

% % Plot mode shapes
% figure('Name','3-DOF Mode Shapes','NumberTitle','off');
% xnodes = 1:3;
% for i=1:3
%     subplot(1,3,i)
%     plot(xnodes, phi3(:,i),'o-','LineWidth',1.5); grid on;
%     title(sprintf('Mode %d, f=%.3f Hz', i, freq3(i)));
%     xlabel('node'); ylabel('Normalized amplitude');
%     ylim([-1.2 1.2])
% end

%% PART B: Discrete chain convergence (varying DOF)
maxN = 20;
modes_to_plot = 3;
freqs_chain = nan(maxN,modes_to_plot);

for N=3:maxN
    Mchain = mc*eye(N);
    % assemble stiffness for N masses connected in series with end springs to ground
    Kc = zeros(N);
    for i=1:N
        if i==1
            Kc(i,i) = kc + kc;
            if N>1, Kc(i,i+1) = -kc; end
        elseif i==N
            Kc(i,i) = kc + kc;
            Kc(i,i-1) = -kc;
        else
            Kc(i,i) = 2*kc;
            Kc(i,i-1) = -kc;
            Kc(i,i+1) = -kc;
        end
    end
    [V,D] = eig(Kc, Mchain);
    [d_sorted, id] = sort(diag(D));
    w = sqrt(d_sorted);
    freqs_chain(N, :) = (w(1:modes_to_plot)/(2*pi))';
end

% plot first three frequencies vs N
%figure('Name','Discrete Chain: freq vs number of masses','NumberTitle','off');
% hold on; grid on;
% Ns = 1:maxN;
% for mode = 1:modes_to_plot
%     %plot(3:maxN, freqs_chain(3:maxN,mode),'o-','LineWidth',1.2);
% end
% xlabel('Number of masses (N)'); ylabel('Frequency (Hz)');
% title('Discrete chain: first 3 natural frequencies vs N');
% legend('Mode 1','Mode 2','Mode 3','Location','NorthEast');
% xlim([3 maxN])
% hold off;

%% PART C: Beam FEM (simply supported) - element matrices & assembly

Ke = @(Le) (E*I/Le^3) * [ 12, 6*Le, -12, 6*Le;
                              6*Le, 4*Le^2, -6*Le, 2*Le^2;
                              -12, -6*Le, 12, -6*Le;
                              6*Le, 2*Le^2, -6*Le, 4*Le^2 ];
Mce = @(Le) (rho*A*Le/420) * [156, 22*Le, 54, -13*Le;
                                 22*Le, 4*Le^2, 13*Le, -3*Le^2;
                                 54, 13*Le, 156, -22*Le;
                                 -13*Le, -3*Le^2, -22*Le, 4*Le^2 ];
Mle = @(Le) (rho*A*Le/2) * eye(4);

% Function to assemble global K and M for given number of elements Ne and mass model
assembleBeam1 = @(Ne, massType) assembleBeamMatrices(Ne, L, Ke, Mce, Mle, massType);

% We'll perform modal analysis for a sample mesh (Ne=10) and then run convergence 2..20
Ne_sample = 10;
[Kg_sample, Mg_cons_sample, nodeCoords] = assembleBeam1(Ne_sample, 'consistent');
[~, Mg_lump_sample, ~] = assembleBeam1(Ne_sample, 'lumped');

% Apply simply-supported BC: nodes 1 and Nnodes -> w=0 (vertical DOF)
% DOF ordering per node: [w_i; theta_i] -> global DOF index = 2*(i-1)+[1 2]
Nnodes = Ne_sample + 1;
totalDOF = 2*Nnodes;
% constrained DOFs (w at nodes 1 and Nnodes)
constrained = [1, 2*(Nnodes-1)+1]; % indices for vertical deflections at end nodes

% reduce matrices
freeDOFs = setdiff(1:totalDOF, constrained);

Krc = Kg_sample(freeDOFs, freeDOFs);
Mrc_cons = Mg_cons_sample(freeDOFs, freeDOFs);
Mrc_lump = Mg_lump_sample(freeDOFs, freeDOFs);

% solve eigenproblem, take first nModes
opts.maxit = 500;
[Vcons, Dcons] = eigs(Krc, Mrc_cons, nModes, 'smallestabs');
[Vlump, Dlump] = eigs(Krc, Mrc_lump, nModes, 'smallestabs');

omega_cons = sqrt(abs(diag(Dcons)));
omega_lump = sqrt(abs(diag(Dlump)));

% Get displacement DOFs (w) from eigenvectors to plot w(x) (extract only translational DOFs)
% Need to reconstruct full vector with zeros at constrained DOFs for plotting
modes_to_plot_show = 3;
figure('Name','Beam mode shapes (simply supported) - Consistent vs Lumped','NumberTitle','off');
xplot = linspace(0,L,200);
for im = 1:modes_to_plot_show
    % expand to full DOF vector
    full_phi_cons = zeros(totalDOF,1);
    full_phi_cons(freeDOFs) = Vcons(:,im);
    w_nodes_cons = full_phi_cons(1:2:totalDOF); % take w DOFs: positions 1,3,5,...
    % same for lumped
    full_phi_lump = zeros(totalDOF,1);
    full_phi_lump(freeDOFs) = Vlump(:,im);
    w_nodes_lump = full_phi_lump(1:2:totalDOF);
    % Interpolate along beam using cubic Hermitian shape functions per element for smooth plotting
    subplot(1,modes_to_plot_show,im)
    hold on; grid on; title(sprintf('Mode %d: f_cons=%.3f Hz, f_lump=%.3f Hz', im, omega_cons(im)/(2*pi), omega_lump(im)/(2*pi)));
    xlabel('x (m)'); ylabel('w (normalized)')
    for e=1:Ne_sample
        % node indexes
        n1 = e; n2 = e+1;
        x1 = (n1-1)*(L/Ne_sample); x2 = (n2-1)*(L/Ne_sample);
        Le = x2 - x1;
        % Hermitian shape functions to evaluate pointwise
        xi = linspace(0,1,20);
        N1 = 1 - 3*xi.^2 + 2*xi.^3;
        N2 = (xi - 2*xi.^2 + xi.^3)*Le;
        N3 = 3*xi.^2 - 2*xi.^3;
        N4 = (-xi.^2 + xi.^3)*Le;
        w_e_cons = N1'*w_nodes_cons(n1) + N2'*full_phi_cons(2*n1) + N3'*w_nodes_cons(n2) + N4'*full_phi_cons(2*n2);
        w_e_lump = N1'*w_nodes_lump(n1) + N2'*full_phi_lump(2*n1) + N3'*w_nodes_lump(n2) + N4'*full_phi_lump(2*n2);
        xp = x1 + xi*(Le);
        plot(xp, w_e_cons / max(abs(w_nodes_cons)),'-','LineWidth',1);
        plot(xp, w_e_lump / max(abs(w_nodes_lump)),'--','LineWidth',1);
    end
    if im==1
        legend('Consistent','Lumped');
    end
    hold off
end

%% =========================
%% PART D: Convergence study for beam (2..20 elements)
%% =========================
Ne_min = 2; Ne_max = 20;
Nevals = Ne_min:Ne_max;
num_elems = numel(Nevals);
modes_cons_freq = nan(num_elems,3);
modes_lump_freq = nan(num_elems,3);
% analytical for simply supported: omega_n = (n^2*pi^2)*sqrt(EI/(rho*A*L^4))
analytical_omega = @(n) (n^2 * pi^2) * sqrt(E*I/(rho*A*L^4));

for ii = 1:num_elems
    Ne = Nevals(ii);
    [Kg, Mg_cons, ~] = assembleBeam1(Ne, 'consistent');
    [~, Mg_lump, ~] = assembleBeam(Ne, 'lumped');
    Nnodes = Ne + 1; totalDOF = 2*Nnodes;
    constrained = [1, 2*(Nnodes-1)+1]; % w at ends = 0
    freeDOFs = setdiff(1:totalDOF, constrained);
    % reduce
    Kred = Kg(freeDOFs, freeDOFs);
    Mcons_red = Mg_cons(freeDOFs, freeDOFs);
    Mlump_red = Mg_lump(freeDOFs, freeDOFs);
    % solve for first 6 modes (safety)
    ksolve = 6;
    [V1,D1] = eigs(Kred, Mcons_red, ksolve, 'smallestabs'); % consistent
    [V2,D2] = eigs(Kred, Mlump_red, ksolve, 'smallestabs'); % lumped
    omega_c = sort(sqrt(abs(diag(D1))));
    omega_l = sort(sqrt(abs(diag(D2))));
    % store first three
    modes_cons_freq(ii,:) = (omega_c(1:3)/(2*pi))';
    modes_lump_freq(ii,:) = (omega_l(1:3)/(2*pi))';
end

% Plot frequency vs number of elements for first 3 modes, compared with analytical
figure('Name','Beam: frequency vs number of elements','NumberTitle','off');
for mode = 1:3
    subplot(3,1,mode)
    hold on; grid on;
    plot(Nevals, modes_cons_freq(:,mode), 'o-', 'LineWidth',1.2)
    plot(Nevals, modes_lump_freq(:,mode), 's--', 'LineWidth',1.2)
    % analytical (constant wrt Ne)
    plot(Nevals, (analytical_omega(mode)/(2*pi))*ones(size(Nevals)), 'k-','LineWidth',1)
    xlabel('Number of elements'); ylabel(sprintf('Mode %d freq (Hz)',mode))
    title(sprintf('Mode %d: convergence (Consistent vs Lumped vs Analytical)', mode));
    legend('Consistent','Lumped','Analytical','Location','Best');
    hold off;
end

% Percentage error vs elements (for first 3 modes)
figure('Name','Beam: % error vs elements','NumberTitle','off');
for mode=1:3
    err_cons = 100*(modes_cons_freq(:,mode) - analytical_omega(mode)/(2*pi)) ./ (analytical_omega(mode)/(2*pi));
    err_lump = 100*(modes_lump_freq(:,mode) - analytical_omega(mode)/(2*pi)) ./ (analytical_omega(mode)/(2*pi));
    subplot(3,1,mode)
    hold on; grid on;
    plot(Nevals, err_cons, 'o-','LineWidth',1.2)
    plot(Nevals, err_lump, 's--','LineWidth',1.2)
    xlabel('Number of elements'); ylabel('% error')
    title(sprintf('Mode %d: percent error vs elements', mode))
    legend('Consistent','Lumped','Location','Best');
    hold off;
end

%% =========================
%% Supporting function definitions
%% =========================
function [Kg, Mg, nodeCoords] = assembleBeamMatrices(Ne, L, KeFunc, McFunc, MlumpedFunc, massType)
    % Assemble global stiffness K and mass M for Euler-Bernoulli beam
    % Ne: number of elements
    % L: total length
    % KeFunc, McFunc, MlumpedFunc: function handles returning element matrices for given Le
    % massType: 'consistent' or 'lumped'
    Le = L/Ne;
    Nnodes = Ne + 1;
    totalDOF = 2*Nnodes;
    Kg = zeros(totalDOF);
    Mg = zeros(totalDOF);
    nodeCoords = linspace(0,L,Nnodes)';
    for e = 1:Ne
        Ke = KeFunc(Le);
        Mc = McFunc(Le);
        Ml = MlumpedFunc(Le);
        % global DOF indices for element e: nodes e and e+1
        n1 = e; n2 = e+1;
        edof = [2*n1-1, 2*n1, 2*n2-1, 2*n2];
        Kg(edof, edof) = Kg(edof, edof) + Ke;
        if strcmpi(massType,'consistent')
            Mg(edof, edof) = Mg(edof, edof) + Mc;
        else
            Mg(edof, edof) = Mg(edof, edof) + Ml;
        end
    end
end

function [Kg, Mg, nodeCoords] = assembleBeam(Ne, L, KeFunc, McFunc, MlumpedFunc, massType)
    % wrapper to match earlier anonymous usage
    [Kg, Mg, nodeCoords] = assembleBeamMatrices(Ne, L, KeFunc, McFunc, MlumpedFunc, massType);
end
