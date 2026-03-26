clear;close all;clc;
A=pi/4 *(0.2*0.2-0.14*0.14);
E=210e9;
I=pi/64 *(0.2^4-0.14^4);
P=10e3;
% Fe=[0,0,0,0,P,0]';

Lv=2.0;
Lh=1.0;
L=Lv+Lh;
Neh=10;
Nev=20;

tdof=3*(Neh+Nev+1);
F=zeros(tdof,1);
F(end-5:end)=Fe;
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

types=[1,2];
thita=[pi/2,0];
for t=types
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
    keU=((A*Le^2)/I)*[1,-1;-1,1];
    keVT =(E*I/Le^3) * [ 12, 6*Le, -12, 6*Le;
                 6*Le, 4*Le^2, -6*Le, 2*Le^2;
                  -12, -6*Le, 12, -6*Le;
                 6*Le, 2*Le^2, -6*Le, 4*Le^2 ];
    Ke=[keU(1,1),0,0,keU(1,2),0,0;
         0,keVT(1,1),keVT(1,2),0,keVT(1,3),keVT(1,4);
         0,keVT(2,1),keVT(2,2),0,keVT(2,3),keVT(2,4);
         keU(2,1),0,0,keU(2,2),0,0;
         0,keVT(3,1),keVT(3,2),0,keVT(3,3),keVT(3,4);
         0,keVT(4,1),keVT(4,2),0,keVT(4,3),keVT(4,4)];
    for e = st:N
            edof = 3*e-2:3*e+3;
            KG(edof, edof) = KG(edof, edof) + T' * Ke * T;
            % MgC(edof, edof) = MgC(edof, edof) + Mce;
            % MgL(edof, edof) = MgL(edof, edof) + Mle;
    end
end

% boundary conditions

fixedDofs = 1:3; 
KG(fixedDofs, :) = [];
KG(:, fixedDofs) = [];
F(fixedDofs)=[];
%fprintf('Global Stiffness matrix');
%disp(KG);
Ug=KG\F;
U=[0;0;0;Ug];
uh=1000*U(3*Nev+1:3:end-2);
uv=1000*U(1:3:3*Nev+1);
vh=U(3*Nev+2:3:end-1);
vv=U(2:3:3*Nev+2);
% fprintf('Displacement in x\n');
% disp(u);
% fprintf('Displacement in y\n');
% disp(v);
% Ul=T'*U(1:6);
% fprintf('Displacement local matrix');
xh=linspace(0,Lh,Neh+1);
xv=linspace(0,0,Nev+1);
yh=linspace(Lv,Lv,Neh+1);
yv=linspace(0,Lv,Nev+1);
xhn=xh+uh';
xvn=xv+uv';
yhn=yh+vh';
yvn=yv+vv';
plot(xh,yh);hold on
plot(xv,yv);
plot(xhn,yhn);
ploy(xvn,yvn);hold off

