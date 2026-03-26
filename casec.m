clear;close all;clc;
N=4;    % number of nodes
K=[100,150,300,400];
Kg=zeros(N);
F=[0,50,0,100]';
NCM=[1,1,2,100;2,2,3,150;3,2,3,300;4,3,4,400];
for i=1:4
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
fe=K(1)*(Um(j)-Um(i));

% Parameters
L = 10;              % Total length of the spring
n_coils = 35;        % Number of coils
amplitude = 0.2;     % Radius of the spring (amplitude of sine wave)
points_per_coil = 500;

% Generate x values before elongation
x1 = linspace(0, L/3, n_coils * points_per_coil);
x21 = linspace(L/3, 2*L/3, n_coils * points_per_coil);
x22 = linspace(L/3, 2*L/3, n_coils * points_per_coil);
x3 = linspace(2*L/3, L, n_coils * points_per_coil);
% Generate x values after elongation
xd1 = linspace(0, L/3+Um(2), n_coils * points_per_coil);
xd21 = linspace(L/3+Um(2), 2*L/3+Um(3), n_coils * points_per_coil);
xd22 = linspace(L/3+Um(2), 2*L/3+Um(3), n_coils * points_per_coil);
xd3 = linspace(2*L/3+Um(3), L+Um(4), n_coils * points_per_coil);



% Generate y values for sine wave
y1 = amplitude * sin(2 * pi * n_coils * x1 / L);
y2 = amplitude * sin(2 * pi * n_coils * x21 / L);
y3 = amplitude * sin(2 * pi * n_coils * x3 / L);



% Plot the spring
figure;hold on;
plot(x1, y1, 'b', 'LineWidth', 2);
plot(x21, y2-1, 'r', 'LineWidth', 2);
plot(x22, y2+1, 'm', 'LineWidth', 2);
plot(x3, y3, 'g', 'LineWidth', 2);
plot(xd1, y1-5, 'b', 'LineWidth', 2);
plot(xd21, y2-6, 'r', 'LineWidth', 2);
plot(xd22, y2-4, 'm', 'LineWidth', 2);
plot(xd3, y3-5, 'g', 'LineWidth', 2);
x2=[L/3,L/3];
x3=[2*L/3,2*L/3];
xd2=[L/3+Um(2),L/3+Um(2)];
xd3=[2*L/3+Um(3),2*L/3+Um(3)];
y=[-1,1];
yd=[-4,-6];
plot(x2,y);
plot(x3,y);
plot(xd2,yd);
plot(xd3,yd);
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