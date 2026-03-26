% FEM Practical 5 - Cantilever Beam Deflection 

clear; clc; close all
E = 200e9;
b = 0.04;
d = 0.08;
I = (b*d^3)/12;
L = 2.0;

P = -1000;
q = -1000;
N=16;% No. of elements
Ne = 1:N;

%% Analytical Solutions

WLp = P*L^3/(3*E*I);

WLu = q*L^4/(8*E*I);
fprintf('Theory Tip dflection for point load\n');
disp(WLp)
fprintf('Theory Tip dflection for UDL\n');
disp(WLu)
resPoint = [];
resUdlLumped = [];
resUdlConsis = [];
Upf=[];
Ulf=[];
Ucf=[];
for Ne = Ne
    Le = L/Ne;                    
    dof = 2*(Ne+1);                
    K = zeros(dof);               
    Fp = zeros(dof,1);      
    FLump = zeros(dof,1);       
    FConsi = zeros(dof,1);    

    ke = (E*I/Le^3) * ...
        [12,   6*Le,  -12,   6*Le;
         6*Le, 4*Le^2,-6*Le, 2*Le^2;
        -12,  -6*Le,   12,   -6*Le;
         6*Le, 2*Le^2,-6*Le, 4*Le^2];

    % Element load vectors
    feLumped = q*Le*[0.5;0;0.5;0];
    feConsis = q*Le*[0.5;Le/12;0.5;-Le/12];

    for e = 1:Ne
        idx = 2*e-1:2*e+2;
        K(idx,idx) = K(idx,idx) + ke;
        FLump(idx) = FLump(idx) + feLumped;
        FConsi(idx) = FConsi(idx) + feConsis;
    end

    Fp(end-1) = P; 

    % Apply boundary conditions (fixed at node 1: w1=0, theta1=0)
    fixedDof = [1,2];
    freeDof = setdiff(1:dof,fixedDof);
    Up = zeros(dof,1);
    ULumped = zeros(dof,1);
    UConsis = zeros(dof,1);

    Up(freeDof) = K(freeDof,freeDof)\Fp(freeDof);
    ULumped(freeDof) = K(freeDof,freeDof)\FLump(freeDof);
    UConsis(freeDof) = K(freeDof,freeDof)\FConsi(freeDof);
     Upf{Ne}={Up};
     Ulf{Ne}={ULumped};
     Ucf{Ne} ={UConsis};
    
    Wlp = Up(end-1);
    Wll = ULumped(end-1);
    Wlc = UConsis(end-1);

   
    resPoint = [resPoint; Ne, Wlp,abs((Wlp-WLp)/WLp*100)];
    resUdlLumped = [resUdlLumped; Ne, Wll,abs((Wll-WLu)/WLu*100)];
    resUdlConsis = [resUdlConsis; Ne, Wlc,abs((Wlc-WLu)/WLu*100)];

    %% Plot deformed shape
    x_nodes = linspace(0,L,Ne+1);
    figure(1); hold on;
    plot(x_nodes,Up(1:2:end),'o-','DisplayName',['Ne=',num2str(Ne)]);
    title('Cantilever with Point Load - Deformed Shape'); xlabel('x (m)'); ylabel('Deflection (m)');
    legend('Location','best'); grid on;

    figure(2); hold on;
    plot(x_nodes,UConsis(1:2:end),'o-','DisplayName',['Ne=',num2str(Ne)]);
    title('Cantilever with UDL (Consistent) - Deformed Shape'); xlabel('x (m)'); ylabel('Deflection (m)');
    legend('Location','best'); grid on;
end
%% Mesh convergence plots


figure(3);
semilogx(resPoint(:,1),resPoint(:,2),'o-','LineWidth',1.5);
hold on;
yline(WLp,'--r','Theory');
xlabel('Number of Elements'); ylabel('Tip Deflection (m)');
title('Mesh Convergence - Point Load'); grid on;



figure(4);
semilogx(resUdlConsis(:,1),resUdlConsis(:,2),'o-','LineWidth',1.5);
hold on;
yline(WLu,'--r','Theory');
xlabel('Number of Elements'); ylabel('Tip Deflection (m)');
title('Mesh Convergence - UDL (Consistent Load)'); grid on;
figure(5);
semilogx(resUdlLumped(:,1),resUdlLumped(:,2),'o-','LineWidth',1.5);
hold on;
yline(Wll,'--r','Theory');
xlabel('Number of Elements'); ylabel('Tip Deflection (m)');
title('Mesh Convergence - UDL (Lumped Load)'); grid on;

%% Display Results
disp('Results: Point Load at Free End');
disp(array2table(resPoint,'VariableNames',{'Ne','TipDeflection','Error_%'}));

disp('Results: UDL (Lumped Load Vector)');
disp(array2table(resUdlLumped,'VariableNames',{'Ne','TipDeflection','Error_%'}));

disp('Results: UDL (Consistent Load Vector)');
disp(array2table(resUdlConsis,'VariableNames',{'Ne','TipDeflection','Error_%'}));


 % t=table(resUdlConsis);
 % writetable(t, 'UDL_Consitance.xlsx');
% t=table(resUdlLumped);
% writetable(t, 'UDL_Lumped.xlsx');
% t = table(resPoint);
% writetable(t, 'Point_Load_Results.xlsx');
