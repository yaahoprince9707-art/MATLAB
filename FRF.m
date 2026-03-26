n=10;
m=1;
k=1;
c=0;
M=m*eye(n);
K=zeros(n);
C=zeros(n);
for row=1:n
    for j=1:row+1
        if row==j
            K(row,j)=2*k;
            C(row,j)=2*c;
          
            if j==n
               break;
            end
        end
        
        if (j==row+1 || j==row-1)
            K(row,j)=-1*k;
            C(row,j)=-1*c;
        end
    end
end
% Frequency range
w = linspace(0, 5,150);
H = zeros(n, length(w));

% Compute frequency response function (FRF)
for i = 1:length(w)
    omega_i = w(i);
    A = -omega_i^2 * M + 1i * omega_i * C + K;
    H(:, i) = inv(A) * eye(n, 1); 
end

% Plot frequency response
% for i=1:n
% figure;
% plot(w, H(i, :), 'LineWidth', 1.5); 
% xlabel('Frequency (rad/s)');
% ylabel('Amplitude');
% title('FRF of mass',{i});
% grid on;
% end
plot(w, H(1, :), 'LineWidth', 1.5); 
xlabel('Frequency (rad/s)');
ylabel('Amplitude');
title('FRF of mass',{1});
grid on;