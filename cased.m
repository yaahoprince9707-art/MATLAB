clear;close all;clc;
N=4;    % number of nodes
K=[100,150,200,400,50,100];
Kg=zeros(N);
F=[0,0,0,100]';
NCM=[1,1,2,100;2,1,3,150;3,1,3,200;4,2,4,400;5,2,4,50;6,3,4,100];
for i=1:6
    r=NCM(i,2:3);
    c=NCM(i,2:3);
    Kg(r,c)=Kg(r,c)+NCM(i,end)*[1,-1;-1,1];
end
Kg;
% After boundary condition
% 1st node is fixed ie. u1=0 is taken as BC
Km=Kg(2:end,2:end);
Fm=F(2:end);  
Um = [0;Km\Fm];
R=Kg*Um-F;  % reaction forces
disp('Reaction forces R:');
disp(R);
i=1;  % 0<i<=N
j=4;  %i<j<=N
Keq=Keqv(Kg,i,F,j);
fe=Keq*(Um(j)-Um(i));
function Ke=Keqv(Kg,ith_node,F,jth)
         Kg(ith_node,:)=[];
         Kg(:,ith_node)=[];
      %Km=Kg(2:end,2:end);
       F(ith_node,:)=[];
        U = [0;Kg\F];
        Ke=(1)/(U(jth));
end


% Parameters
L = 10;              % Total length of the spring
n_coils = 35;        % Number of coils
amplitude = 0.2;     % Radius of the spring (amplitude of sine wave)
points_per_coil = 500;

% Generate x values
x1 =  linspace(0, L/2, n_coils * points_per_coil);
x2 =  linspace(0, L/2, n_coils * points_per_coil);
x3 =  linspace(0, L/2, n_coils * points_per_coil);
x4 = linspace(L/2, L, n_coils * points_per_coil);
x5 = linspace(L/2, L, n_coils * points_per_coil);
x6 = linspace(L/2, L, n_coils * points_per_coil);

xd1 =  linspace(0, L/2+Um(2), n_coils * points_per_coil);
xd2 =  linspace(0, L/2+Um(3), n_coils * points_per_coil);
xd6 =  linspace(L/2+Um(3), L+Um(4), n_coils * points_per_coil);
xd4 = linspace(L/2+Um(2), L+Um(4), n_coils * points_per_coil);

% Generate y values for sine wave
y1 = amplitude * sin(2 * pi * n_coils * x1 / L);
y2 = amplitude * sin(2 * pi * n_coils * x2 / L);
y3 = amplitude * sin(2 * pi * n_coils * x3 / L);



% Plot the spring
figure;hold on;
plot(x1, y1, 'b', 'LineWidth', 2);
plot(x2, y2-2, 'r', 'LineWidth', 2);
plot(x3, y3-4, 'g', 'LineWidth', 2);
plot(x4, y1+1, 'm', 'LineWidth', 2);
plot(x5, y2-1, 'c', 'LineWidth', 2);
plot(x6, y3-3, 'color','black ','LineWidth', 2);
plot(xd1, y1-7, 'b', 'LineWidth', 2);
plot(xd2, y2-9, 'r', 'LineWidth', 2);
plot(xd2, y3-11, 'g', 'LineWidth', 2);
plot(xd4, y1-6, 'm', 'LineWidth', 2);
plot(xd4, y2-8, 'c', 'LineWidth', 2);
plot(xd6, y3-10, 'color','black ','LineWidth', 2);

x2=[L/2,L/2];
x3=[L/2,L/2];
x4=[L,L];

xld2=[L/2+Um(2),L/2+Um(2)];
xld3=[L/2+Um(3),L/2+Um(3)];
xld4=[L+Um(4),L+Um(4)];
xl=[0,0];
yl=[-11,0];
plot(xl,yl,'Color','black');
yld2=[-8,-6];
yld3=[-11,-9];
yld4=[-10,-6];
yl2=[-1,1];
yl3=[-4,-2];
yl4=[-3,1];
plot(x2,yl2);
plot(x3,yl3);
plot(x4,yl4);
plot(xld2,yld2);
plot(xld3,yld3);
plot(xld4,yld4);

% plot(xd2,yd);
% plot(xd3,yd);
axis equal;
xlabel('Length');
ylabel('Displacement');
title('Spring Representation');
% xlim([min(x)-1, max(x)+1]);
grid on;
hold off;





