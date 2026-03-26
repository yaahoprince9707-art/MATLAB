clear;close all;clc;
b=[0.02,0.010];
E=[205e9,70e9];
Mz=500;
fx=2000;
fy=3000;
%Fe=[0,0,0,0,P,0]';

Lv=2.0;
Lh=1.0;
L=Lv+Lh;
Neh=10;
Nev=20;
Nnodes=Nev+Neh+1;
tdof=3*(Neh+Nev+1);
F=zeros(tdof,1);
KG=zeros(tdof);
%R=[cos(thita),-sin(thita),0;sin(thita),cos(thita),0;0,0,1];
Tr = @(c, s) ...
    [ c, s, 0, 0, 0, 0;
     -s, c, 0, 0, 0, 0;
      0, 0, 1, 0, 0, 0;
      0, 0, 0, c, s, 0;
      0, 0, 0, -s, c, 0;
      0, 0, 0, 0, 0, 1 ];

% T=blkdiag(R,R);
% K = T' * Ke * T;

types=[1,2];% 1 for vertical & 2 for horizontal
thita=[pi/2,0];
for t=types
    A=b(t)^2;
    I=(b(t)^4)/12;
    c = cos(thita(t));
    s = sin(thita(t));
    T = Tr(c, s);
    st=1;
    Le=Lv/Nev;
    N=Nev;
    if t==2
        st=Nev+1;
        N=Neh+Nev;
        Le=Lh/Neh;
    end
    keA = (E(t)*A/Le)*[1 -1; -1 1];
    keB = (E(t)*I/Le^3)*[ 12 6*Le -12 6*Le;
                        6*Le 4*Le^2 -6*Le 2*Le^2;
                       -12 -6*Le 12 -6*Le;
                        6*Le 2*Le^2 -6*Le 4*Le^2];
    Ke = zeros(6);
    Ke([1 4],[1 4]) = keA;
    Ke([2 3 5 6],[2 3 5 6]) = keB;
    for e = st:N
            edof = 3*e-2:3*e+3;
            KG(edof, edof) = KG(edof, edof) + T' * Ke * T;
            if e==Nev+1
                 % F(edof(1))=fx;
                 %  F(edof(2))=fy;
                 F(edof(3))=-Mz;
            end
            % MgC(edof, edof) = MgC(edof, edof) + Mce;
            % MgL(edof, edof) = MgL(edof, edof) + Mle;
    end
end
% boundary conditions

fixed = [1:3,edof(4:6)]; 
free = setdiff(1:tdof,fixed);
U = zeros(tdof,1);
U(free) = KG(free,free)\F(free);
% fprintf('Global Stiffness matrix');
% disp(KG);

% Ua=[0;0;0;U;0;0;0];
% disp(Ua);
X = zeros(Nnodes,1); Y = zeros(Nnodes,1);
for i=1:(Nev+1)
    X(i)=0; Y(i)=(i-1)*(Lv/Nev);
end
for i=1:Neh
    idx=Nev+1+i;
    X(idx)=-i*(Lh/Neh); Y(idx)=Lv;
end

U3 = reshape(U,3,[])';
u=U3(:,1); v=U3(:,2);

scale = 1;
Xd = X + scale*u;
Yd = Y + scale*v;
%plots
hold on; %axis equal; grid on; box on;
set(gca,'FontSize',11);
plot(X, Y, 'r-','LineWidth',1.2,'MarkerFaceColor','k','DisplayName','Undeformed');
plot(Xd, Yd, 'b-.','LineWidth',1.5,'MarkerFaceColor','r','DisplayName','Deformed');
legend('Location','best');
xlabel('X (m)'); ylabel('Y (m)');
title(sprintf('Deformed L-Frame (scale = %d×)',scale));