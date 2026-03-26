clear all;
clc;
It=1:20;
Iy=function1(It);
num_int=length(It);
VM=[];
for i=1:num_int
    VM=[VM;It(i).^(0:num_int-1)];
end
W=inv(VM)*Iy';
polyNomial=@(x) x.^(0:num_int-1)*W;
x=linspace(0,20,20);
y=function1(x);
z=zeros(1,length(x));
for i=1:length(x)
    z(i)=polyNomial(x(i));
end
figure
grid on
hold on;
plot(x,y);
%plot(It,Iy,'o');
%plot(x,z);
hold off;


