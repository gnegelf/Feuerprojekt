function [ A,b_U,b_L,c,solInhom,solBasis,Aext,b_Uext,b_Lext,Aext2,b_Uext2,b_Lext2] = constrContamination(i1,i2,i3, xn,tn,dx,dt,N,u,contPos,step,C,Tmax,Schwell,p1,p2)
%This function uses diriclet boundary conditions, which might be
%inappropriate
paramsControlled=p1;
paramsInhom=p2;
xn1=xn+1;
tn1=tn+1;
i4=(xn+1)^2*(tn+1);
N1=N+1;

contBasis=zeros(i4,N); %Anzahl der Controller
global usetime;
global time;
gridPoints=zeros((xn+1)^2,2);
for yy=0:xn
    for xx=0:xn
       gridPoints(xx*(xn+1)+yy+1,1)=xx/xn;
       gridPoints(xx*(xn+1)+yy+1,2)=yy/xn;
    end
end
gridPoints=[gridPoints;contPos];
gridPoints=round(gridPoints,8);
solBasis=zeros(i4,i2);
solBasis2=zeros(size(contPos,1)*(tn+1),i2);
zet=zeros(i4,N);
zet2=zeros(size(contPos,1)*(tn+1),N);
gridParams=struct('LL',[0,0],'UR',[1,1],'h_max',1/20);

for i=1:N
   %u=@(x) ofem.matrixarray(-10*exp(-50*dot(x-contPos(i,:)',x-contPos(i,:)',1)));
   T_xt=heatRobin(gridPoints,0:dt:time,@(x) u(x,contPos(i,:)'),@(x) ofem.matrixarray(0+0*x(1,1,:)),paramsControlled,gridParams);
   %Eliminate boundary conditions!
   store=reshape(T_xt(1:(xn+1)^2,:),(xn+1)^2*size(T_xt,2),1);
   zet2(:,i)=reshape(T_xt((xn+1)^2+1:size(T_xt,1),:),size(contPos,1)*size(T_xt,2),1);
   zet(:,i)=store;
end

solBasis(:,1:N)=zet;
solBasis2(:,1:N)=zet2;
for t=2:tn
     %solStore(:,(t-1)*N+1:t*N)=[zeros((t-1)*xn1^2,N);zet(xn1^2+1:(tn1-t+2)*xn1^2,:)];
     solBasis(:,(t-1)*N+1:t*N)=[zeros((t-1)*xn1^2,N);zet(1:(tn1-t+1)*xn1^2,:)];
     solBasis2(:,(t-1)*N+1:t*N)=[zeros((t-1)*size(contPos,1),N);zet2(1:(tn1-t+1)*size(contPos,1),:)];
end

solInhomRaw=heatRobin(gridPoints,0:dt:time,@(x) 0*ofem.matrixarray(x(1,1,:)),p1.initialT,paramsInhom,gridParams);
solInhom=reshape(solInhomRaw(1:(xn+1)^2,:),(xn+1)^2*size(solInhomRaw,2),1);
solInhom2=reshape(solInhomRaw(1+(xn+1)^2:size(solInhomRaw,1),:),size(contPos,1)*size(solInhomRaw,2),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Networkconditions
A=sparse(tn*N+1,tn*N+N);

for t=1:tn
    for i=1:N
        A((t-1)*N+i,tn*N+i)=-1;
        A((t-1)*N+i,(t-1)*N+i)=1;
    end
end
for i=1:N
   A(tn*N+1,tn*N+i)=1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate b_L and b_U
%b_U=zeros(i3+i2+i2,1);
%b_U=zeros(i3+i2+i2+count,1);
b_U=zeros(tn*N+1,1);
b_L=-100*ones(tn*N+1,1);
b_L(tn*N+1)=0;
b_U(tn*N+1)=5;


Aext= [solBasis sparse(i1,N)];
b_Lext= -solInhom;
b_Uext= 1000000.1*ones(i1,1);

Aext2= [solBasis2 sparse(size(contPos,1)*(tn+1),N)];
b_Lext2= -solInhom2;
b_Uext2= 1000000.1*ones(size(contPos,1)*(tn+1),1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%construct c
n=size(A,2);
c=zeros(n,1);
factors=zeros(tn1*(xn+1)^2,1);

for t=1:tn1
    if (t==1)
       fact=0.5;
    else
        fact=1;
    end
    for j=1:xn+1
        if (j==1 || j==xn+1)
           factors((t-1)*(xn+1)^2+j)=0.25*fact; 
        else
            factors((t-1)*(xn+1)^2+j)=0.5*fact;
        end
    end
end
c(1:i2)=solBasis'*factors;
end