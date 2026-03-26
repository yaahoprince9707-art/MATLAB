clear;clc;
n=32;          % No. of elements 
m=1:n+1;      % No. of Nodes or mesh points
L=1000;
P=50000;
F=zeros(n+1,1);
F(end)=P;
d1=40;d2=20;E=210000;
x=linspace(0,L,n+1);
Ax=@(x) (pi/4)*((d1+(d2-d1)*(x/L).^2)).^2;

Kga=zeros(n+1);
Kgm=zeros(n+1);
dx=@(x) ((pi/4)*((d1+(d2-d1)*(x/L).^2)).^2).^(-1);


Analyt_def=(P/E)*integral(dx,0,L);
fprintf('The analytical deflection is: %.6f\n', Analyt_def);
Aea=zeros(n,1);  % Calculate avg area
Aem=zeros(n,1);  % Calculate mid-point area
for i=1:n
    Le =x(i+1)-x(i); % Length of each element
    Aea(i)=(1/Le)*integral(Ax,x(i),x(i+1)); % Calculate avg area for each segment
    Aem(i) =Ax((x(i)+x(i+1))/2);            % Calculate mid-point area for each segment
end

c0l=1:n;c01=c0l(:);
c02=c01;
c03=c02+1;
c04a=Aea;
c04m=Aem;
NCAa=[c01,c02,c03,c04a];
NCAm=[c01,c02,c03,c04m];
for i=1:n
    ra=NCAa(i,2:3);
    ca=NCAa(i,2:3);
    rm=NCAa(i,2:3);
    cm=NCAa(i,2:3);
    Kgm(ra,ca)=Kgm(rm,cm)+Ke(Aem(i),x(i+1)-x(i));
    Kga(ra,ca)=Kga(ra,ca)+Ke(Aea(i),x(i+1)-x(i));
end
disp(Kga);
disp(Kgm);
% After boundary condition
% 1st node is fixed ie. u1=0 is taken as BC
Kma=Kga(2:end,2:end);
Kmm=Kgm(2:end,2:end);
Fm=F(2:end);  
Uma = [0;Kma\Fm];
Umm = [0;Kmm\Fm];
% Calculate the displacement vector
disp('Calculated displacements as per exact average:');
disp(Uma);
disp('Calculated displacements as per mid-point area:');
disp(Umm);
% Error Calculation ----------
 err_ea=(Analyt_def-Uma)/Analyt_def *100;
 err_mp=(Analyt_def-Umm)/Analyt_def *100;
disp(err_ea);
disp(err_mp);


 figure;
    plot(x, zeros(size(x)), 'k--', 'LineWidth', 1.5,'DisplayName','original'); hold on;
    plot(x, Uma, 'b-o', 'DisplayName', 'Exact Avg'); 
    plot(x, Umm, 'r-*', 'DisplayName', 'Mid-point');
    xlabel('Position along bar (m)');
    ylabel('Displacement (m)');
    title(['Deformed Shape for nelem = ', num2str(n)]);
    legend; grid on;
% cross section of the bar
x_plot = linspace(0, L, 100);
d_plot = d1 + (d2 - d1)*(x_plot/L).^2;
A_plot = pi/4 * d_plot.^2;
figure;
plot(x_plot, A_plot*1e6, 'm', 'LineWidth', 2);
xlabel('Position along bar (mm)');
ylabel('Cross-sectional Area (mm^2)');
title('Cross-section Profile');
grid on;

figure;
hold on
semilogy(m,err_ea,'r-o');
semilogy(m,err_mp,'b--');
title('Nelem Vs Error');
xlabel('Nelem');
ylabel('error');
legend('error exact avg','error mid-point')
hold off;
grid on;

function Keq=Ke(Ae,Le)
E=210000;
Keq = ((E * Ae) / Le)*[1,-1;-1,1];
end
