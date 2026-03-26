clear; close all; clc;
m=5;
k=100;
Ne=10;
modestp = 3;
freqs = nan(Ne,modestp);
for N=3:Ne
    Mc = m*eye(N);
    Kc = zeros(N);
    for i=1:N
        if i==1
            Kc(i,i) = k + k;
            if N>1, Kc(i,i+1) = -k; end
        elseif i==N
            Kc(i,i) = k + k;
            Kc(i,i-1) = -k;
        else
            Kc(i,i) = 2*k;
            Kc(i,i-1) = -k;
            Kc(i,i+1) = -k;
        end
    end
    [U,val] = eig(Kc, Mc);
    [wSorted, id] = sort(diag(val));
    w = sqrt(wSorted);
    freqs(N, :) = (w(1:modestp)/(2*pi))';
end

% plot first three frequencies vs N
styl={'-';'-.';'--'};
mark={'o','s','d'};
colours=lines(3);
figure;
hold on; grid on;

for mode = 1:modestp
    plot(3:Ne, freqs(3:Ne,mode),'Color', colours(mode,:), ...
        'LineStyle', styl{mode}, ... 
        'Marker', mark{mode}, ...  
        'LineWidth', 1.5, ...
        'MarkerSize', 6,'LineWidth',1.2);
end
xlabel('Number of masses (N)'); ylabel('Frequency (Hz)');
title('first 3 natural frequencies vs N');
legend('Mode 1','Mode 2','Mode 3','Location','NorthEast');
xlim([3 Ne])
hold off;

