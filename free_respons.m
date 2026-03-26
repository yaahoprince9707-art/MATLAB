n=5;
m=1;
k=1;
M=m*eye(n);
K=zeros(n);
for row=1:n
    for j=1:row+1
        if row==j
            K(row,j)=2*k;
          
            if j==n
               break;
            end
        end
        
        if (j==row+1 || j==row-1)
            K(row,j)=-1*k;
        end
    end
end
K(n,n)=1;
[U , val]=eig(K,M);
W=diag(sqrt(val));
for row=1:n
    U(:,row)=U(:,row)/U(1,row);
end
disp(U)
disp(W)
% we are having M, K, W, and mode shaps till line-26
% initial conditions matrix
X0=rand(n,1);
V0=rand(n,1);
Y0=[X0;V0];
% formation of initial 
IC=zeros(n*2);
m2=1; %starting no. of row of IC matrix
for row=1:2
    m1=1;
    for j=1:2*n
        if rem(j,2)==0
            m1=m1-1;
        end
        IC(m2:row*n,j)=(-1)^((j-1)*(row+1))*U(:,m1)*(W(m1)*1i)^(row-1);
        m1=m1+1;
    end
    m2=n+1;
end
disp(IC)