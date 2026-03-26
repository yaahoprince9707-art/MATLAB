% cross section of the bar
x_plot = linspace(0, L, 100);
y_plot=x_plot;
d_plot = d1 + (d2 - d1)*(x_plot/L).^2;
A_plot = pi/4 * d_plot.^2;
B_plot=-A_plot;
figure;
plot(x_plot, A_plot*1e6, 'm', 'LineWidth', 2);
plot(y_plot, B_plot*1e6, 'g', 'LineWidth', 2);
xlabel('Position along bar (mm)');
ylabel('Cross-sectional Area (mm^2)');
title('Cross-section Profile');
grid on;
