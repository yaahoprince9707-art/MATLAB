x = 0:0.1:10;
y = sin(x);
data = [x' y'];
writematrix(data, 'plot_data.xlsx');
