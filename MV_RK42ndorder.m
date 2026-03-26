k=1000;
m=10;
c=50;
h=0.01;
t=0:h:30;

y=zeros(2,length(t));
 y(:,1)=[0.5 0];
for n=1:length(t)-1
    k1=h*fx(y(:,n),k,m,c,t(n));
    k2=h*fx(y(:,n)+k1/2,k,m,c,t(n));
    k3=h*fx(y(:,n)+k2/2,k,m,c,t(n));
    k4=h*fx(y(:,n)+k3,k,m,c,t(n));
    kf=(1/6)*(k1+2*k2+2*k3+k4);
    y(:,n+1)=y(:,n)+kf;
    %disp(y(:,n+1));
end

hold on
plot(t,y(1,:));
plot(t,y(2,:));
hold off
% plot(y(1,:),y(2,:));
 title('phase protrat')
 grid on

function dydt=fx(y,k,m,c,t)
ft=[0;(1/m)*1.5*exp(pi*t*1i)];
dydt=[0 1 ;(-k/m) (-c/m)]*y+ft;
end