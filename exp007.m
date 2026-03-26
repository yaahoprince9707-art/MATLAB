clc; clear; close all;

do = 0.20; di = 0.14;
A = pi/4*(do^2 - di^2);
E = 210e9;
I = pi/64*(do^4 - di^4);

Lv = 2.5; Lh = 0.5;
Nev = 20;  Neh = 5;
Nnodes = Nev + Neh + 1;
tdof = 3 * Nnodes;
F = zeros(tdof,1);
P = -1e4;
F(end-1) = P;
Tr = @(c,s) [ c, -s, 0, 0, 0, 0;
              s,  c, 0, 0, 0, 0;
              0,  0, 1, 0, 0, 0;
              0,  0, 0, c, -s, 0;
              0,  0, 0, s,  c, 0;
              0,  0, 0, 0,  0, 1 ];

KG = zeros(tdof);
types  = [1 2];
theta  = [pi/2 0];

for t = types
    c = cos(theta(t)); s = sin(theta(t)); T = Tr(c,s);
    if t==1, st=1; N=Nev; Le=Lv/Nev;
    else,    st=Nev+1; N=Nev+Neh; Le=Lh/Neh; end
    
    keA = (E*A/Le)*[1 -1; -1 1];
    keB = (E*I/Le^3)*[ 12 6*Le -12 6*Le;
                        6*Le 4*Le^2 -6*Le 2*Le^2;
                       -12 -6*Le 12 -6*Le;
                        6*Le 2*Le^2 -6*Le 4*Le^2];
    Ke = zeros(6);
    Ke([1 4],[1 4]) = keA;
    Ke([2 3 5 6],[2 3 5 6]) = keB;
    
    for e = st:N
        edof = 3*(e-1)+1 : 3*(e+1);
        KG(edof,edof) = KG(edof,edof) + T'*Ke*T;
    end
end
% BCs
fixed = 1:3;
free = setdiff(1:tdof,fixed);

U = zeros(tdof,1);
U(free) = KG(free,free)\F(free);
%plots
X = zeros(Nnodes,1); Y = zeros(Nnodes,1);
for i=1:(Nev+1)
    X(i)=0; Y(i)=(i-1)*(Lv/Nev);
end
for i=1:Neh
    idx=Nev+1+i;
    X(idx)=i*(Lh/Neh); Y(idx)=Lv;
end

U3 = reshape(U,3,[])';
u=U3(:,1); v=U3(:,2);

scale = 100;
Xd = X + scale*u;
Yd = Y + scale*v;

hold on; %axis equal; grid on; box on;
set(gca,'FontSize',11);
plot(X, Y, '-b','LineWidth',1.2,'MarkerFaceColor','k','DisplayName','Undeformed');
plot(Xd, Yd, '-.r','LineWidth',1.5,'MarkerFaceColor','r','DisplayName','Deformed');
legend('Location','best');
xlabel('X (m)'); ylabel('Y (m)');
title(sprintf('Deformed L-Frame (scale = %d×)',scale));
