function [ A,b_U,b_L,c,solInhom,solBasis,Aext,b_Uext,b_Lext,Aext2,b_Uext2,b_Lext2] = constrMIPEliminationRobinVertexFE(i1,i2,i3, xn,tn,dx,dt,N,u,contPos,step,C,Tmax,Schwell,p1,p2)
%This function uses diriclet boundary conditions, which might be
%inappropriate
paramsControlled=p1;
paramsInhom=p2;
xn1=xn+1;
tn1=tn+1;
i4=(xn+1)^2*(tn+1);
N1=N+1;
contBasis=zeros(i4,N); %Anzahl der Controller
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

T_ints=zeros(N,tn);
for i=1:N
   %u=@(x) ofem.matrixarray(-10*exp(-50*dot(x-contPos(i,:)',x-contPos(i,:)',1)));
   [T_xt,T_ints(i,:)]=heatRobinNew(gridPoints,0:dt:time,@(x) u(x,contPos(i,:)'),@(x) ofem.matrixarray(0+0*x(1,1,:)),paramsControlled,gridParams);
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


[TM,tempArcInhom]=calcMaxCondEliminationVertex(i1,i2,i3, contPos,xn,tn,N,C,solBasis,solInhom,Tmax);
s2=2*i2+i3;
condM=sparse(i3+i2+i2,s2);
adderC=i2;

%%%%%%%%%%%%%%%%%%Kapazit�tenbedingungen
for t=1:tn+1
    for i=1:N+1
        for j=1:N+1
            condM((t-1)*N1^2+(i-1)*N1+j,i2+(t-1)*N1^2+(i-1)*N1+j)=1;
            
            if (i>1) && (j>1)%the input from source is independent of the heat flow
                if t>1
                    condM((t-1)*N1^2+(i-1)*N1+j,i2+i3+(t-2)*N+j-1)=-C(i,j,1);
                end
                %condM((t-1)*N1^2+(i-1)*N1+j,i2+2*i3+(t-1)*N1^2+(i-1)*N1+j)=-C(i,j,1);    %Kapazit�tenbedingung MIP
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%Flowconditions
for t=2:tn+1
    for i=1:N
        condM(i3+(t-2)*N+i,(t-2)*N+i)=1; %set coefficient of variable corresponding to wateroutput at controller i at time t-1
        for j=1:N+1
            condM(i3+(t-2)*N+i,i2+(t-1)*N1^2+(i)*N1+j)=1;
            
            if t>step*C(i+1,j,2)
                condM(i3+(t-2)*N+i,i2+(t-step*C(i+1,j,2)-1)*N1^2+(j-1)*N1+i+1)=-1;
            end
        end
    end   
end

for t=1:tn
    for i=1:N
        condM(i3+i2+(t-1)*N+i,(t-1)*N+i)=1;
        condM(i3+i2+(t-1)*N+i,i3+i2+(t-1)*N+i)=-1;
    end
end
%count=0;
% for t=2:tn+1
%     for i=1:N
%         for j=1:N
%             if C(j,i,1) > 0
%                 count=count+1;
%                 condM(i3+2*i2+count,i2+(t-1)*(N+1)^2+j*(N+1)+i+1)=1;
%                 condM(i3+2*i2+count,i3+i2+(t-2)*N+i)=-1;
%             end
%             
%         end
%     end
% end

A=condM;
A=[A;TM];
V=sparse(1,2*i2+i3);

%start flow = zero condition
for i=1:N
    V(1,i)=1;
    for j=1:N+1
        V(1,i2+(i)*N1+j)=1;
    end
end
A=[A;V];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate b_L and b_U
%b_U=zeros(i3+i2+i2,1);
%b_U=zeros(i3+i2+i2+count,1);
b_U=zeros(i3+i2+i2,1);
for t=1:tn+1
    for j=1:N+1
       b_U((t-1)*(N+1)^2+j)=C(1,j,1);
       b_U((t-1)*(N+1)^2+(j-1)*(N+1)+1)=C(j,1,1);
    end
end
b_L=zeros(size(b_U));
b_L(1:i3)=-1000*ones(i3,1);
b_L(i3+1:i3+i2)=-10000*ones(i2,1);
%b_L(i3+i2+1:i3+i2+i2+count)=-1*ones(i2+count,1);
b_L(i3+i2+1:i3+i2+i2)=-1*ones(i2,1);
b_U=[b_U;-tempArcInhom+Tmax+Schwell];
b_L=[b_L;-tempArcInhom];

%start flow zeros rhs
b_U=[b_U;0];
b_L=[b_L;0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Stateconstraints
 %A=[A; [solBasis sparse(i1,3*i3)]];
 %b_L=[b_L; -solInhom];
 %b_U=[b_U; 1000000*ones(i1,1)];
%Aext= [[solBasis;solBasis2] sparse(i1+size(contPos,1)*(tn+1),i3+i2)];
%b_Lext= -[solInhom;solInhom2];
%b_Uext= 1000000.1*ones(i1+size(contPos,1)*(tn+1),1);

Aext= [solBasis sparse(i1,i3+i2)];
b_Lext= -solInhom;
b_Uext= 1000000.1*ones(i1,1);

Aext2= [solBasis2 sparse(size(contPos,1)*(tn+1),i3+i2)];
b_Lext2= -solInhom2;
b_Uext2= 1000000.1*ones(size(contPos,1)*(tn+1),1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%construct c
n=size(A,2);
c=zeros(n,1);
factors=dt*dx*dx*ones(tn1*(xn+1)^2,1);

for t=1:tn1
    if (t==1 || t==tn1)
       fact=0.5; 
    else
        fact=1;
    end
    for i=[1,xn+1]
        for j=1:xn+1
            if (j==1 || j==xn+1)
               factors((t-1)*(xn+1)^2+(i-1)*(xn+1)+j)=dt*dx*dx*0.25*fact; 
            else
                factors((t-1)*(xn+1)^2+(i-1)*(xn+1)+j)=dt*dx*dx*0.5*fact;
                factors((t-1)*(xn+1)^2+(j-1)*(xn+1)+i)=dt*dx*dx*0.5*fact;
            end
        end
    end
end
c(1:i2)=solBasis'*factors;
c2=zeros(i2,1);
for t=1:tn
   for i=1:N
       c2((t-1)*N+i)=sum(T_ints(i,1:tn+1-t));
   end
end
c(1:i2)=c2;
end
