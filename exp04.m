clear;clc;
N=3;
n=2*N;
L=1;
P=1000;
F=[0,0,0,0,0,-P]';
ke=1000*[1,-1;-1,1];
KG=zeros(n);
NCA=[1,1,2;2,2,3;3,3,1];
cs=T(pi/3,L,1,N);
Tr=[cs(2),cs(1),0,0;0,0,cs(2),cs(1)];
Kg1=Tr'*ke*Tr;
cs=T(pi/3,L,2,N);
Tr=[cs(2),cs(1),0,0;0,0,cs(2),cs(1)];
Kg2=Tr'*ke*Tr;
cs=T(pi/3,L,3,N);
Tr=[cs(2),cs(1),0,0;0,0,cs(2),cs(1)];
Kg3=Tr'*ke*Tr;
Kg={Kg1,Kg2,Kg3};
for i=1:3
    n=[NCA(i,2),NCA(i,3)];
    r=[2*n(1)-1,2*n(1),2*n(2)-1,2*n(2)];
    c=[2*n(1)-1,2*n(1),2*n(2)-1,2*n(2)];
    KG(r,c)=KG(r,c)+Kg{i};
end

% After boundary condition
% 1st node is fixed ie. u1=0 is taken as BC
KGa=KG;
Fa=F;
KG(1:4,:)=[];
KG(:,1:4)=[];
F=F(5:end);
disp(KG);
U=KG\F;
fprintf('Nodal Displacements\n');

U=[0;0;0;0;U];
disp(U);
% Calculate the reaction forces at the supports
reactionForces = KGa * U - Fa;
disp('Reaction Forces at Supports');
disp(reactionForces);
% Plotting the undeformed shape
%fig = figure('Color','w','Name','Truss-Deformation');
figure
hold on;
title('Deformed and Undeformed Shape of the Truss Elements');
xlabel('X-axis');
ylabel('Y-axis');

% Define original node positions
originalNodes = [0, 0; 1, 0; cos(pi/3), sin(pi/3); 0, 0; 0, 0]; % Adjust based on actual node positions
deformedNodes = originalNodes + [U(1), U(2); U(3), U(4); U(5), U(6); 0, 0; 0, 0]; % Apply displacements

% Plot undeformed shape
plot(originalNodes(:,1), originalNodes(:,2), 'b--.', 'DisplayName', 'Undeformed Shape');

% Plot deformed shape
plot(deformedNodes(:,1), deformedNodes(:,2), 'r-o', 'DisplayName', 'Deformed Shape');

legend show;
grid on;
axis equal;
hold off;

function Trans=T(thita,L,e,N)
r1=e+1;
if e==N
    r1=1;
end
node=[0,0;1,0;cos(thita),sin(thita)];

s=(node(r1,2)-node(e,2))/L;
c=(node(r1,1)-node(e,1))/L;
Trans=[s,c];
end
