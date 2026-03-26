%% p-type refinement
clear;clc;
n=50;
p=10;q=1:n;
error=zeros(1,length(q));
x=linspace(0,2,1000);
y=function2(x);
fx=@(x) function2(x);
Aa=integral(fx,0,2);
fprintf('Actual area under the curve: %.8f\n', Aa);
   %Method of list squares to fit polynomials
for i=1:n
% Discretize the x-axis
A=vanderM(x,i);
c=(inv(A'*A))*A'*y'; 
% Compute area using trapz
area=calc_Area(c,x(1),x(end));
error(i) = abs(((Aa - abs(area))/Aa)*100);  % populating errors
fprintf('Approximate area under the curve using %d degree polynomial is: %.6f\n', i, area);
end
% area=zeros(1,length(q));
% dA=0;i=2;
% A=vanderM(x,1);
% c=(inv(A'*A))*A'*y';
% area(1) = calc_Area(c, x(1), x(end)); % Initialize the first area calculation
% while dA<0.01
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
plot(x,y);
plot(x,z);
legend(sprintf('y = f(x)'), sprintf('polynomial of degree %d', p));
hold off;
figure;
semilogy(q,error);
title('degree of polynomials Vs % error');
xlabel('p');
ylabel('error');
