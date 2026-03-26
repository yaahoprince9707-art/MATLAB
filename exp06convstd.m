clc;clear;close all;

E = 200e9;
b = 0.04;
d = 0.08;
I = (b*d^3)/12;
L = 2.0;
Nemin=2;
Nemax=20;
Ne=Nemin:Nemax;
r=7850;
a=b*d;
freqc=nan(numel(Ne),3);
freql=nan(numel(Ne),3);
Wt = @(n) (n^2 * pi^2) * sqrt(E*I/(r*a*L^4));
for e = 1:numel(Ne)
    Le = L / Ne(e);
    A=(r*a*Le)/420;
    Mle =(r*A*Le/2) * eye(4);
    Mce=A*[156,22*Le,54,-13*Le; 
        22*Le, 4*Le*Le ,13*Le ,-3*Le*Le;
        54,13*Le,156,-22*Le;
        -13*Le,-3*Le*Le,-22*Le,4*Le*Le];
    Ke =(E*I/Le^3) * [ 12, 6*Le, -12, 6*Le;
                 6*Le, 4*Le^2, -6*Le, 2*Le^2;
                  -12, -6*Le, 12, -6*Le;
                 6*Le, 2*Le^2, -6*Le, 4*Le^2 ];
    N = Ne(e) + 1;
    tDOF = 2*N;
    Kg = zeros(tDOF);
    MgC = zeros(tDOF);
    MgL = zeros(tDOF);
    for i=1:Ne(e)
         n1 = i; n2 = i+1;
         edof = [2*n1-1, 2*n1, 2*n2-1, 2*n2];
         Kg(edof, edof) = Kg(edof, edof) + Ke;
         MgC(edof, edof) = MgC(edof, edof) + Mce;
         MgL(edof, edof) = MgL(edof, edof) + Mle;
    end
    fxdof = [1, 2*Ne(e)+1];
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
    Wl2=diag(valL);
    Wl2=real(Wl2);
    Wl = sqrt(Wl2);
    freqc(e,:) = (Wc(1:3)/(2*pi))';
    freql(e,:) = (Wl(1:3)/(2*pi))';
end

for mode = 1:3
    figure
    hold on; grid on;
    plot(Ne, freqc(:,mode), 'o-', 'LineWidth',1.2)
    plot(Ne, freql(:,mode), 's--', 'LineWidth',1.2)
    % analytical (constant wrt Ne)
    plot(Ne, (Wt(mode)/(2*pi))*ones(size(Ne)), 'k-','LineWidth',1)
    xlabel('Number of elements'); ylabel(sprintf('Mode %d freq (Hz)',mode))
    title(sprintf('Mode %d: convergence (Consistent vs Lumped vs Analytical)', mode));
    legend('Consistent','Lumped','Analytical','Location','Best');
    hold off;
end


for mode=1:3
    errc = 100*(freqc(:,mode) - Wt(mode)/(2*pi)) ./ (Wt(mode)/(2*pi));
    errl = 100*(freql(:,mode) - Wt(mode)/(2*pi)) ./ (Wt(mode)/(2*pi));
    % disp('%age error in natural freqc.');
    % disp(array2table(errc,'VariableNames',{'Error_%'}));
    % disp('%age error in %d natural freql.');
    % disp(array2table(errl,'VariableNames',{'Error_%'}));
    figure
    hold on; grid on;
    plot(Ne, errc, 'o-','LineWidth',1.2)
    plot(Ne, errl, 's--','LineWidth',1.2)
    xlabel('Number of elements'); ylabel('% error')
    title(sprintf('Mode %d: percent error vs elements', mode))
    legend('Consistent','Lumped','Location','Best');
    hold off;
end