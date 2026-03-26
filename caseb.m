clear;close all;clc;
N=2;    % number of nodes
K=[80,50,200];
Kg=zeros(N);
F=[0,80]';
NCM=[1,1,2,80;2,1,2,50;3,1,2,200];
for i=1:3
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
i=2;  % 0<i<=N
j=2;  %i<j<=N
Keq=Keqv(Kg,i,F,j);
fe=K(3)*(Um(j)-Um(1));


% Parameters
L = 10;              % Total length of the spring
n_coils = 35;        % Number of coils
amplitude = 0.2;     % Radius of the spring (amplitude of sine wave)
points_per_coil = 500;

% Generate x values before elongation
x1 = linspace(0, L, n_coils * points_per_coil);
x2 = linspace(0, L, n_coils * points_per_coil);
x3 = linspace(0, L, n_coils * points_per_coil);
% Generate x values after elongation
xd1 = linspace(0, L+Um(2), n_coils * points_per_coil);
xd2 = linspace(0, L+Um(2), n_coils * points_per_coil);
xd3 = linspace(0, L+Um(2), n_coils * points_per_coil);



% Generate y values for sine wave
y1 = amplitude * sin(2 * pi * n_coils * x1 / L);
y2 = amplitude * sin(2 * pi * n_coils * x2 / L);
y3 = amplitude * sin(2 * pi * n_coils * x3 / L);
yd1 = amplitude * sin(2 * pi * n_coils * xd1 / L);
yd2 = amplitude * sin(2 * pi * n_coils * xd2 / L);
yd3 = amplitude * sin(2 * pi * n_coils * xd3 / L);



% Plot the spring
figure;hold on;
plot(x1, y1, 'b', 'LineWidth', 2);
plot(x2, y2-1, 'r', 'LineWidth', 2);
plot(x3, y3-2, 'g', 'LineWidth', 2);
plot(xd1, yd1-4, 'b', 'LineWidth', 2);
plot(xd2, yd2-5, 'r', 'LineWidth', 2);
plot(xd3, yd3-6, 'g', 'LineWidth', 2);
xl1=[L,L];
yl1=[-2,0];
plot(xl1,yl1);
xld1=[L+Um(2),L+Um(2)];
yld1=[-6,-4];
plot(xld1,yld1);
axis equal;
xlabel('Length');
ylabel('Displacement');
title('Spring Representation');
% xlim([min(x)-1, max(x)+1]);
grid on;
hold off;


function Ke=Keqv(Kg,ith_node,F,jth)
         Kg(ith_node,:)=[];
         Kg(:,ith_node)=[];
      %Km=Kg(2:end,2:end);
       F(ith_node,:)=[];
        U = [0;Kg\F];
        Ke=(1)/(U(jth));

end