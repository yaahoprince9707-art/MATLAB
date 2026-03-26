n=5;
m=1;
k=1;
c=1;
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
K(n,n)=1;
C(n,n)=1;
h=0.01;
t=0:0.01:100;
y=zeros(2*n,length(t));
y(:,1)=rand(2*n,1);
for m=1:length(t)-1
    k1=h*fx(y(:,m),K,M,C,t(m),n);
    k2=h*fx(y(:,m)+k1/2,K,M,C,t(m),n);
    k3=h*fx(y(:,m)+k2/2,K,M,C,t(m),n);
    k4=h*fx(y(:,m)+k3,K,M,C,t(m),n);
    kf=(1/6)*(k1+2*k2+2*k3+k4);
    y(:,m+1)=y(:,m)+kf;
    %disp(y(:,n+1));
end


for i=1:2*n
    figure
    xt=y(i,:);
    plot(t,xt);
    title('respons',{i});
    grid on
end
% plot(t,y(1,:),t,y(2,:));
%  grid on
function dydt=fx(y,K,M,C,t,n)
A00=zeros(n);A11=eye(n);
An1=zeros(n,1);
f=zeros(n,1);
f(n,1)=1.5*exp(pi*t*1i);
ft=[An1;(-inv(M))*f];
dydt=[A00 A11;-inv(M)*K -inv(M)*C]*y+ft;
end
