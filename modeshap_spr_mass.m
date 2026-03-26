clear; close all; clc;

n=3;
m=5;
M=m*eye(n);
k=100;
K=zeros(n);
for i=1:n
    for j=1:i+1
        if i==j
            K(i,j)=2*k;
            if j==n
               break;
            end
        end
        if (j==i+1 || j==i-1)
            K(i,j)=-1*k;
        end
    end
end
%K(n,n)=1;
[U , val]=eig(K,M);
[wsort, id] = sort(diag(val));
w = sqrt(wsort);
freqs = w/(2*pi)';
fprintf('Natural frequencies (Hz)\n')
disp(freqs)
for i=1:n
    U(:,i)=U(:,i)/U(1,i);
end

for i=1:n
    figure
    y=[0;U(:,i);0];
    plot(y);
    title('mode shap',{i});
    grid on
end
