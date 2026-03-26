clc;
clear all;
k1=4000;k2=2000;k3=6000;
m1=1;m2=2;
M=[m1 0; 0 m2];
K=[k1+k2 -k2; -k2 k2+k3];
[modesp ,frq]=eig(K,M);
A00=zeros(2);A11=eye(2);
CM=[A00 A11;-inv(M)*K A00];
max_frq=max(sqrt(diag(frq))/2*pi);
dt=1/(max_frq*20);
t=0:dt:3000*dt;
testode=@(t,y)(CM*y);

y0=[0.01 0 0 0];
[time ,yt]=ode23(testode,t,y0);
subplot(4, 1,1)
plot(t,yt(:,1),'r');
subplot(4, 1,2)
plot(t,yt(:,2),'b')
subplot(4, 1,3)
plot(t,yt(:,3),'g')
subplot(4, 1,4)
plot(t,yt(:,4),'o')
%xlim(-0.02,0.02)
%ylim(-0.02,0.02)

