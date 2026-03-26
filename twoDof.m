clc;
clear all;
k=4;
k1=k;k2=2*k;k3=k;
m1=1;m2=1;
M=[m1 0; 0 m2];
K=[k1+k2 -k2; -k2 k2+k3];
k11=k1+k2;
k22=k2+k3;
k12=-k2;
x10=0.05;x20=1;v10=0;v20=0;

y0=[x10;x20;v10;v20];% initial condition matrix
p=[1 -((k11/m1)+(k22/m2)) (k11*k22-k12*k12)/m1*m2];% natural frequency eqn
r=roots(p);%w2^2 and w1^2
w1w2=sqrt(r);
% time step
max_frq=max(w1w2);
dt=1/(max_frq*20);
t=0:dt:3000*dt;
u1=[1;(-k12)/(k22-r(2,1)*m2)];%mode shap 1
u2=[1;(-k12)/(k22-r(1,1)*m2)];% mode shap2
iw1u1=w1w2(2)*u1*1i;
iw2u2=w1w2(1)*u2*1i;
IC=[u1 u1 u2 u2;iw1u1 -iw1u1 iw2u2 -iw2u2];
C14=inv(IC)*y0;
disp(C14);
y1=zeros(1,length(t));
y2=zeros(1,length(t));
v1=zeros(1,length(t));
v2=zeros(1,length(t));
function yt=x1t(t,IC,C14)
yt=exp(t*IC(3,:))*C14;
end
function yt=x2t(t,IC,C14)
yt=IC(2,:).*exp(t*IC(3,:))*C14;
end
function yt=v1t(t,IC,C14)
yt=IC(3,:).*exp(t*IC(3,:))*C14;
end
function yt=v2t(t,IC,C14)
yt=IC(4,:).*exp(t*IC(3,:))*C14;
end
for n=1:length(t)
    y1(1,n)=x1t(t(1,n),IC,C14);
    y2(1,n)=x2t(t(1,n),IC,C14);
    v1(1,n)=v1t(t(1,n),IC,C14);
    v2(1,n)=v2t(t(1,n),IC,C14);
end
subplot(4, 1,1)

plot(t,y1,'r');
title('x1')
subplot(4, 1,2)

plot(t,y2,'b')
title('x2')
subplot(4, 1,3)

plot(t,v1,'g')
title('v1')
subplot(4, 1,4)

plot(t,v2,'cyan')
title('v2')