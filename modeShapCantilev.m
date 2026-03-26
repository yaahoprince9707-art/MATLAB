
clear; close all; clc;

L = 2.0;              
ne = 40;              
nn = ne + 1;         
le = L / ne;          

E = 210e9;            
rho = 7850;            
b = 0.02;            
h = 0.02;             
A = b * h;             
I = (b * h^3) / 12;    
q = 1000;              


ndof = 2 * nn;         
K = zeros(ndof);
M = zeros(ndof);
F = zeros(ndof,1);    


for e = 1:ne
    n1 = e;
    n2 = e + 1;
    dofs = [2*n1-1, 2*n1, 2*n2-1, 2*n2]; 
    ke = E*I / le^3 * [  12,   6*le,  -12,   6*le;
                         6*le, 4*le^2, -6*le, 2*le^2;
                        -12,  -6*le,   12,  -6*le;
                         6*le, 2*le^2, -6*le, 4*le^2 ];
                     
    me = rho*A*le/420 * [156, 22*le, 54, -13*le;
                         22*le, 4*le^2, 13*le, -3*le^2;
                         54, 13*le, 156, -22*le;
                        -13*le, -3*le^2, -22*le, 4*le^2 ];
    fe = q * le / 12 * [6; le; 6; -le]; 
    K(dofs,dofs) = K(dofs,dofs) + ke;
    M(dofs,dofs) = M(dofs,dofs) + me;
    F(dofs) = F(dofs) + fe;
end

fixedDOF = [1, 2];
allDOF = 1:ndof;
freeDOF = setdiff(allDOF, fixedDOF);

Kff = K(freeDOF, freeDOF);
Ff = F(freeDOF);
uf = Kff \ Ff;

u = zeros(ndof,1);
u(freeDOF) = uf;

wStatic = u(1:2:end); 
Krr = K(freeDOF, freeDOF);
Mrr = M(freeDOF, freeDOF);

nModes = 3;
opts.maxit = 500;
[phiR, D] = eigs(Krr, Mrr, nModes, 'SM', opts);
omega2 = diag(D);
omega = sqrt(real(omega2));
freq = omega / (2*pi);

Phi = zeros(ndof, nModes);
for i = 1:nModes
    Phi(freeDOF, i) = phiR(:, i);
end


xnodes = linspace(0, L, nn)';
modesW = Phi(1:2:end, :);


for i = 1:nModes
    modesW(:, i) = modesW(:, i) / max(abs(modesW(:, i)));
    if mean(modesW(:, i)) < 0
        modesW(:, i) = -modesW(:, i);
    end
end


nm_to_show = min(nModes, 6);
for i = 1:nm_to_show
    
    figure
    plot(xnodes, modesW(:,i), '-o','LineWidth',1.2,'MarkerSize',4);
    hold on;
    scaleStatic = 5 * max(abs(modesW(:,i))) / max(abs(wStatic)+eps);
    plot(xnodes, -scaleStatic*wStatic, '--','LineWidth',1.1);
    hold off;
    title(sprintf('Mode %d — f = %.3f Hz', i, freq(i)));
    xlabel('x (m)');
    ylabel('Normalized displacement');
    legend('Mode shape (FEM)','Scaled static deflection','Location','best');
    grid on;
end
sgtitle('Cantilever Beam — FEM Mode Shapes and Static Deflection (UDL)');
fprintf('Natural frequencies (Hz) for cantilever beam:\n');
for i = 1:nm_to_show
    fprintf(' Mode %d: %.6f Hz\n', i, freq(i));
end
