function [ A,b_U,b_L,c,solInhom,solBasis,Aext,b_Uext,b_Lext] = constrMIPEliminationRobin(i1,i2,i3, R,xn,tn,dx,dt,N,u,contPos,a1,a2,step,C,g,f,Tmax,Schwell,hL,hR,hU,hO,D)
%This function uses diriclet boundary conditions, which might be
%inappropriate
xn1=xn+3;
i4=xn1^2*(tn+1);
N1=N+1;
tn1=tn+1;
%A=sparse(xn1^2*(tn+1),xn1^2*(tn+1));
%A(1:xn1^2,1:xn1^2)=speye(xn1^2,xn1^2);
%R=@(x,y) D(x,y)*dt/dx^2;
R=D*dt/dx^2;
%discretization of Laplace
difM=sparse(xn1^2,2*xn1^2);
%difM=full(difM);

for i=2:xn+2
    for j=2:xn+2
        %nDiag=-0.5*R((i-2)*dx,(j-2)*dx);
        %rx=dt/dx*a1((i-2)*dx,(j-2)*dx);
        %ry=dt/dx*a2((i-2)*dx,(j-2)*dx);
        nDiag=-0.5*R;
        rx=dt/dx*a1;
        ry=dt/dx*a2;
        %difM((i-1)*xn1+j,(i-1)*xn1+j+xn1^2  )=1+2*R((i-2)*dx,(j-2)*dx);
        difM((i-1)*xn1+j,(i-1)*xn1+j+xn1^2  )=1+2*R;
        difM((i-1)*xn1+j,(i-2)*xn1+j+xn1^2  )=nDiag; 
        difM((i-1)*xn1+j,(i-1)*xn1+j-1+xn1^2)=nDiag;
        difM((i-1)*xn1+j,(i-1)*xn1+j+1+xn1^2)=nDiag;
        difM((i-1)*xn1+j, i   *xn1+j+xn1^2  )=nDiag;
        %difM((i-1)*xn1+j,(i-1)*xn1+j  )=-1+2*R((i-2)*dx,(j-2)*dx);
        difM((i-1)*xn1+j,(i-1)*xn1+j  )=-1+2*R;
        if rx>0
            difM((i-1)*xn1+j,(i-2)*xn1+j )=nDiag-rx;
            difM((i-1)*xn1+j, i   *xn1+j)=nDiag;
            difM((i-1)*xn1+j,(i-1)*xn1+j  )=difM((i-1)*xn1+j,(i-1)*xn1+j)+rx;
        else
            difM((i-1)*xn1+j,(i-2)*xn1+j )=nDiag;
            difM((i-1)*xn1+j, i   *xn1+j)=nDiag+rx;
            difM((i-1)*xn1+j,(i-1)*xn1+j  )=difM((i-1)*xn1+j,(i-1)*xn1+j)-rx;
        end
        if ry>0
            difM((i-1)*xn1+j,(i-1)*xn1+j-1)=nDiag;
            difM((i-1)*xn1+j,(i-1)*xn1+j+1)=nDiag-ry;
            difM((i-1)*xn1+j,(i-1)*xn1+j  )=difM((i-1)*xn1+j,(i-1)*xn1+j)+ry;
        else
            difM((i-1)*xn1+j,(i-1)*xn1+j-1)=nDiag+ry;
            difM((i-1)*xn1+j,(i-1)*xn1+j+1)=nDiag;
            difM((i-1)*xn1+j,(i-1)*xn1+j  )=difM((i-1)*xn1+j,(i-1)*xn1+j)-ry;
        end
    end
end
for j=2:xn+2
    difM(j,j+xn1^2  )=1/(2*dx);
    difM(j,2*xn1+j+xn1^2  )=-1/(2*dx);
    difM(j,xn1+j+xn1^2  )=hL;
    difM((xn1-1)*xn1+j,(xn1-1)*xn1+j+xn1^2  )=+1/(2*dx);
    difM((xn1-1)*xn1+j,(xn1-3)*xn1+j+xn1^2  )=-1/(2*dx);
    difM((xn1-1)*xn1+j,(xn1-2)*xn1+j+xn1^2  )=hR;
end
for i=2:xn+2
    difM((i-1)*xn1+1,(i-1)*xn1+1+xn1^2  )=1/(2*dx);
    difM((i-1)*xn1+1,(i-1)*xn1+3+xn1^2  )=1/-(2*dx);
    difM((i-1)*xn1+1,(i-1)*xn1+2+xn1^2  )=hU;
    difM((i-1)*xn1+xn1,(i-1)*xn1+xn1+xn1^2  )=1/(2*dx);    
    difM((i-1)*xn1+xn1,(i-1)*xn1+xn1-2+xn1^2  )=1/-(2*dx);
    difM((i-1)*xn1+xn1,(i-1)*xn1+xn1-1+xn1^2  )=hO;
end
difM(1,1+xn1^2 )=1;
difM((xn1-1)*xn1+1,(xn1-1)*xn1+1+xn1^2  )=1;
difM(xn1,xn1+xn1^2  )=1;
difM((xn1-1)*xn1+xn1,(xn1-1)*xn1+xn1+xn1^2  )=1;
A=sparse(xn1^2*(tn+1),xn1^2*(tn+1));
A(1:xn1^2,1:xn1^2)=speye(xn1^2,xn1^2);
Aadd=kron(speye(tn),difM(:,1:xn1^2));
Aadd2=kron(speye(tn),difM(:,xn1^2+1:xn1^2*2));
A(xn1^2+1:xn1^2*(tn+1),1:xn1^2*(tn))=Aadd;
A(xn1^2+1:xn1^2*(tn+1),xn1^2+1:xn1^2*(tn+1))=A(xn1^2+1:xn1^2*(tn+1),xn1^2+1:xn1^2*(tn+1))+Aadd2;
contBasis=zeros(i4,N); %Anzahl der Controller
for t=2
    for i=2:xn+2
        for j=2:xn+2
            for k=1:N
                if abs(u((i-2)*dx,(j-2)*dx,contPos(k,1),contPos(k,2)))>0.0001
                    contBasis((t-1)*xn1^2+(i-1)*xn1+j,(t-2)*N+k)=u((i-2)*dx,(j-2)*dx,contPos(k,1),contPos(k,2));
                    %A((t-1)*xn1^2+(i-1)*xn1+j,tn1*xn1^2+(t-2)*N+k)=u((j-1)*dx,(i-1)*dx,contPos(k,1),contPos(k,2));%vorher: u((i-1)*dx,(j-1)*dx,contPos(k,1),contPos(k,2));
                end
            end
        end
    end
end
solStore=zeros(i4,i2);
zet=A\contBasis;
solStore(:,1:N)=zet;
for t=2:tn
     %solStore(:,(t-1)*N+1:t*N)=[zeros((t-1)*xn1^2,N);zet(xn1^2+1:(tn1-t+2)*xn1^2,:)];
     solStore(:,(t-1)*N+1:t*N)=[zeros((t-1)*xn1^2,N);zet(1:(tn1-t+1)*xn1^2,:)];
end

b=zeros(i4,1);
for t=1:tn+1
    if t==1
        for k=2:xn+2
            for j=2:xn+2
                b((t-1)*(xn1)^2+(j-1)*(xn1)+k)=g(dx*(j-2),dx*(k-2));
            end
        end
    else
        for k=1:xn1
            for j=1:xn1  
                if j==1 || j==xn1 || k==1 || k==xn1 
                    b((t-1)*(xn1)^2+(j-1)*(xn1)+k)=(k==1)*hL+(j==1)*hU+(k==xn1)*hR+(j==xn1)*hO;
                else

                end
            end
        end
    end
    b((t-1)*(xn1)^2+1)=0;
    b((t-1)*(xn1)^2+xn1)=0;
    b((t-1)*(xn1)^2+(xn1-1)*(xn1)+1)=0;
    b((t-1)*(xn1)^2+(xn1-1)*(xn1)+xn1)=0;
end


solInhomStore=A\b;
solBasis=zeros(i1,i2);
solInhom=zeros(i1,1);
for t=1:tn+1
    for x=1:xn+1
        for y=1:xn+1
            solBasis((t-1)*(xn+1)^2+(x-1)*(xn+1)+y,:)=solStore((t-1)*(xn+3)^2+(x)*(xn+3)+y+1,:);
            solInhom((t-1)*(xn+1)^2+(x-1)*(xn+1)+y)=solInhomStore((t-1)*(xn+3)^2+(x)*(xn+3)+y+1);
        end
    end
end
xn1=xn+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Networkconditions


[TM,tempArcInhom]=calcMaxCondElimination(i1,i2,i3, contPos,xn,tn,N,C,solBasis,solInhom );
s2=i2+i3*3;
condM=sparse(i3+i2,s2);
adderC=i2;

%%%%%%%%%%%%%%%%%%Kapazit�tenbedingungen
for t=1:tn+1
    for i=1:N+1
        for j=1:N+1
            condM((t-1)*N1^2+(i-1)*N1+j,adderC+(t-1)*N1^2+(i-1)*N1+j)=1;
            if (i>1) && (j>1)%the input from source is independent of the heat flow
                condM((t-1)*N1^2+(i-1)*N1+j,i2+2*i3+(t-1)*N1^2+(i-1)*N1+j)=-C(i,j,1);    %Kapazit�tenbedingung MIP
                for cc=0:C(i,j,2)-1 %additional flow that might still be on the arc from earlier timeframes
                    if t>C(i,j,2)*step
                        for pp=0:step-1
                            condM((t-1)*N1^2+(i-1)*N1+j,adderC+(t-1-pp-cc*step)*N1^2+(i-1)*N1+j)=1;
                        end
                    end
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%Flowconditions
for t=2:tn+1
    for i=1:N
        condM(i3+(t-2)*N+i,(t-2)*N+i)=1; %set coefficient of variable corresponding to wateroutput at controller i at time t-1
        for j=1:N+1
            condM(i3+(t-2)*N+i,adderC+(t-1)*N1^2+(i)*N1+j)=1;
            if t>step*C(i+1,j,2)
                condM(i3+(t-2)*N+i,adderC+(t-step*C(i+1,j,2)-1)*N1^2+(j-1)*N1+i+1)=-1;
            end
        end
    end   
end
A=condM;
A=[A;TM];
IntMat=sparse(N^2*tn1,i2+i3*3);
for t=1:tn1
   for i=2:N1
       for j=2:N1
           if C(i,j,1)>0
               if t>1
                    IntMat((t-1)*N^2+N*(i-2)+j-1,i2+2*i3+(t-1)*N1^2+N1*(i-1)+j)=Tmax;
                    IntMat((t-1)*N^2+N*(i-2)+j-1,i2+  i3+(t-1)*N1^2+N1*(i-1)+j)=1;
               else
                    IntMat((t-1)*N^2+N*(i-2)+j-1,i2+2*i3+(t-1)*N1^2+N1*(i-1)+j)=1;
               end
           end

       end
   end
end
A=[A;IntMat];
V=zeros(1,s2);
%start flow = zero condition
for i=1:N
    V(1,i)=1;
    for j=1:N+1
        V(1,adderC+(i)*N1+j)=1;
    end
end
A=[A;V];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate b_L and b_U
s1=size(TM,1);

b_U=zeros(i3+i2,1);
b_L=b_U;
b_L(1:i3)=-1000000.1*ones(i3,1);
for t=1:tn+1
    for i=1:N+1
       b_U((t-1)*(N+1)^2+(i-1)*(N+1)+1)=C(i,1,1);
       b_U((t-1)*(N+1)^2+i)=C(1,i,1);
    end
end
b_U=[b_U;-tempArcInhom];
b_L=[b_L;-1000000.1*ones(s1+(tn+1)*(N)^2,1)];
b_U=[b_U;zeros(N^2,1)];
b_U=[b_U;(Tmax+Schwell)*ones((tn)*(N)^2,1)];
b_U=[b_U;0];
b_L=[b_L;0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Stateconstraints
 %A=[A; [solBasis sparse(i1,3*i3)]];
 %b_L=[b_L; -solInhom];
 %b_U=[b_U; 1000000*ones(i1,1)];
Aext= [solBasis sparse(i1,3*i3)];
b_Lext= -solInhom;
b_Uext= 1000000.1*ones(i1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%construct c
n=size(A,2);
c=zeros(n,1);
c(1:i2)=solBasis'*ones(tn1*(xn+1)^2,1);
c(n-i3+1:n)=-0.00001*ones(i3,1);

end