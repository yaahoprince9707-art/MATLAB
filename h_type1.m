%%--------------------------------- h-type refinement----------------------

clc;
clear;close all;
clc;
a = 0;            % Start of interval
b = 2;            % End of interval
nh = 6;           % Number of trapezoids
m=1:20;
errorh=zeros(1,length(m));

xh = linspace(a, b, 1000);
yh = function1(xh);
Aah=trapz(xh,yh);
fprintf('Actual area under the curve: %.4f\n', Aah);
for i=1:length(m)
% Discretize the x-axis
x_t = linspace(a, b, i+1);
y_t = function1(x_t);
% Compute area using trapz
areah = trapz(x_t, y_t);
errorh(i) = abs(((Aah - areah)/Aah)*100);  % populating errors
fprintf('Approximate area under the curve using %d trapezoids: %.4f\n', i, areah);
end
x_t = linspace(a, b, nh+1);
y_t = function1(x_t);
% Plot the curve
figure;
hold on;
plot(xh, yh, 'LineWidth', 2);  % Original curve
title('Curve y = f(x) with Trapezoidal Approximation');
xlabel('x');
ylabel('f(x)');

%Overlay trapezoids
for i = 1:nh
    % Coordinates of trapezoid
    x_trap = [x_t(i), x_t(i), x_t(i+1), x_t(i+1)];
    y_trap = [0, y_t(i), y_t(i+1), 0];
    fill(x_trap, y_trap, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'r');
end
legend(sprintf('y = f(x)'), sprintf('%d trapezium', nh));
hold off;


% ---------------------------------------------P-TYPE REFINEMENT --------------------------------

n=20;
p=6;q=1:n;
error=zeros(1,length(q));
x=linspace(0,2,20);
y=function1(x);
fx=@(x) function1(x);
Aa=integral(fx,0,2);
fprintf('Actual area under the curve: %.8f\n', Aa);
   % Method of list squares to fit polynomials
for i=1:n
% Discretize the x-axis
A=vanderM(x,i);
c=(inv(A'*A))*A'*y'; 
% Compute area using trapz
area=calc_Area(c,x(1),x(end));
error(i) = abs(((Aa - area)/Aa))*100;  % populating errors
fprintf('Approximate area under the curve using %d degree polynomial is: %.6f\n', i, area);
end
% area=zeros(1,length(q));
% dA=0;i=2;
% A=vanderM(x,1);
% c=(inv(A'*A))*A'*y';
% area(1) = calc_Area(c, x(1), x(end)); % Initialize the first area calculation
% while dA<0.001
%     A = vanderM(x, i);
%     c = (inv(A' * A)) * A' * y';
%     area(i) = calc_Area(c, x(1), x(end));
%     dA=((area(i)-area(i-1))/area(i-1)) *100;
%     i=i+1;
% end
% fprintf('The polynomial for the best fit is having the degree of %d and having area =%.4f:', i-1,area(i-1));


A=vanderM(x,p);
c=(inv(A'*A))*A'*y';

polyNomial=@(x) x.^(0:p)*c;
z=zeros(1,length(x));
for i=1:length(x)
    z(i)=polyNomial(x(i));
end
figure
grid on
hold on;
plot(x,y,'r-o');
plot(x,z);
legend(sprintf('y = f(x)'), sprintf('polynomial of degree %d', p));
hold off;
figure;
subplot(2,1,1)
semilogy(m,errorh);
title('P-TYPE Vs H-TYPE CONVERGENCE');
ylabel('error');
xlabel('n');
subplot(2,1,2)
semilogy(q,error);

xlabel('p');
ylabel('error');

exportgraphics(gcf, 'mySubplots.pdf'); % Saves as PDF

% title('degree of polynomials Vs % error');

