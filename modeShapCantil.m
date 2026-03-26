clc;clear;close all;

E = 200e9;
b = 0.04;
d = 0.08;
I = (b*d^3)/12;
L = 2.0;
Ne=64;
Le = L / (Ne);

roh=7850;
a=b*d;
A=(roh*a*Le)/420;
Mle =(roh*A*Le/2) * eye(4);
Mce=A*[156,22*Le,54,-13*Le; 
    22*Le, 4*Le*Le ,13*Le ,-3*Le*Le;
    54,13*Le,156,-22*Le;
    -13*Le,-3*Le*Le,-22*Le,4*Le*Le];
Ke =(E*I/Le^3) * [ 12, 6*Le, -12, 6*Le;
                 6*Le, 4*Le^2, -6*Le, 2*Le^2;
                  -12, -6*Le, 12, -6*Le;
                 6*Le, 2*Le^2, -6*Le, 4*Le^2 ];
N = Ne + 1;
tDOF = 2*N;
Kg = zeros(tDOF);
MgC = zeros(tDOF);
MgL = zeros(tDOF);
for e = 1:Ne
    n1 = e; n2 = e+1;
    edof = [2*n1-1, 2*n1, 2*n2-1, 2*n2];
    Kg(edof, edof) = Kg(edof, edof) + Ke;
    MgC(edof, edof) = MgC(edof, edof) + Mce;
    MgL(edof, edof) = MgL(edof, edof) + Mle;
end

fxdof = [1, 2];
% modyfied matrices
frDOFs = setdiff(1:tDOF, fxdof);
Kgm = Kg(frDOFs, frDOFs);
MgCm = MgC(frDOFs, frDOFs);
MgLm = MgL(frDOFs, frDOFs);

[Qc , valC]=eig(Kgm,MgCm);
[Ql , valL]=eig(Kgm,MgLm);
Wc2 = diag(valC);
Wc2 = real(Wc2);
Wc = sqrt(Wc2);
freqc = Wc / (2*pi);
Wl2=diag(valL);
Wl2=real(Wl2);
Wl = sqrt(Wl2);
freql = Wl / (2*pi);
Nm=3;
qc = zeros(tDOF, Nm);
ql = zeros(tDOF, Nm);
for i = 1:Nm
    qc(frDOFs,i) = real(Qc(:,i));
    ql(frDOFs,i) = real(Ql(:,i));
end

% extract transverse DOFs for plotting
xn = linspace(0, L, N)';
mwc= qc(1:2:end, :);
mwl= ql(1:2:end, :);

for i = 1:Nm
    mwc(:,i) = mwc(:,i) / max(abs(mwc(:,i)));
    mwl(:,i) = mwl(:,i) / max(abs(mwl(:,i)));
    if mean(mwc(:,i)) < 0
       mwc(:,i) = -mwc(:,i);
    end
    if mean(mwl(:,i)) < 0
       mwl(:,i) = -mwl(:,i);
    end
end


for i = 1:Nm
    % subplot(3,2,i);
    figure
    hold on
    plot(xn, mwc(:,i), 'LineWidth',1.2,'MarkerSize',4);
    plot(xn, mwl(:,i), 'rs','LineWidth',1.2,'MarkerSize',4);
    hold off
    title(sprintf('Mode %d — f = %.3f Hz', i, freqc(i)));
    xlabel('x (m)');
    ylabel('Normalized transverse displacement');
    legend('Mode shape cons','Mode shape lump','Location','best');
    grid on;
end
