function [ A,b_U,b_L,c,solInhom,solBasis ] = constrMIPElimination(i1,i2,i3, R,xn,tn,dx,dt,N,u,contPos,a1,a2,r1,r2,step,C,g,f,Tmax,Schwell)
%This function uses diriclet boundary conditions, which might be
%inappropriate
A=sparse(1,i2+i3*3); %+state constraints
xn1=xn+1;
tn1=tn+1;
N1=N+1;
difM=speye(xn1^2);
difMbig=sparse(i1,i1);
Diag=1+4*R;
nDiag=-R;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PDE calculations
for i=1:xn+1
    for j=1:xn+1
        if (i==1)||(j==1)||(i==xn+1)||(j==xn+1)
            difM((i-1)*(xn+1)+j,(i-1)*(xn+1)+j)=1;
        else
            difM((i-1)*xn1+j,(i-1)*xn1+j  )=Diag;
            difM((i-1)*xn1+j,(i-2)*xn1+j  )=nDiag; 
            difM((i-1)*xn1+j,(i-1)*xn1+j-1)=nDiag;
            difM((i-1)*xn1+j,(i-1)*xn1+j+1)=nDiag;
            difM((i-1)*xn1+j, i   *xn1+j  )=nDiag;
        end  
    end
end
%controlM=[controlM;zeros(N,N)];
%difMext=[difM controlM];
contBasis=zeros(i1,tn*N); %Anzahl der Controller
for t=2:tn+1
    difMbig((t-1)*(xn1^2)+1:t*(xn1^2),(t-1)*(xn1^2)+1:t*xn1^2)=difM;
    for i=2:xn
        for j=2:xn
            for k=1:N
                if u((j-1)*dx,(i-1)*dx,contPos(k,1),contPos(k,2))<0.001
                    contBasis((t-1)*xn1^2+(i-1)*xn1+j,(t-2)*N+k)=u((j-1)*dx,(i-1)*dx,contPos(k,1),contPos(k,2));
                    %A((t-1)*xn1^2+(i-1)*xn1+j,tn1*xn1^2+(t-2)*N+k)=u((j-1)*dx,(i-1)*dx,contPos(k,1),contPos(k,2));%vorher: u((i-1)*dx,(j-1)*dx,contPos(k,1),contPos(k,2));
                end
            end
        end
    end
end
difMbig(1:xn1^2,1:xn1^2)=speye(xn1^2,xn1^2);

%upwind Method
for t=2:tn+1
    for i=1:xn+1
        for j=1:xn+1
            ind =(t-1)*((xn+1)^2)+(i-1)*(xn+1)+j;
            ind1=(t-2)*((xn+1)^2)+(i-1)*(xn+1)+j;
            if ((i==1)||(j==1)||(i==xn+1)||(j==xn+1))
            else
                if a1>0
                    difMbig(ind,ind1)=difMbig(ind,ind1)-1-r1;
                    difMbig(ind,ind1+1)=difMbig(ind,ind1+1)+r1;
                else
                    difMbig(ind,ind1)=difMbig(ind,ind1)-1+r1;
                    difMbig(ind,ind1-1)=difMbig(ind,ind1-1)-r1;
                end
                if a2>0
                    difMbig(ind,ind1)=difMbig(ind,ind1)-r2;
                    difMbig(ind,ind1+xn+1)=difMbig(ind,ind1+xn+1)+r2;
                else
                    difMbig(ind,ind1)=difMbig(ind,ind1)+r2;
                    difMbig(ind,ind1-xn-1)=difMbig(ind,ind1-xn-1)-r2;
                end
            end   
        end
    end
end
solBasis=difMbig\contBasis;
b=zeros(i1,1);
for t=1:tn+1
    for j=1:xn+1
        for k=1:xn+1
            if t==1
                b((t-1)*(xn+1)^2+(j-1)*(xn+1)+k)=g(dx*(j-1),dx*(k-1));
            else
                if j==1 || k==1 || j==xn+1 || k==xn+1
                    b((t-1)*(xn+1)^2+(j-1)*(xn+1)+k)=f(dx*(j-1),dx*(k-1),(t-1)*dt);
                else
                    
                end
            end
        end
    end
end
solInhom=difMbig\b;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Networkconditions


[TM,tempArcInhom]=calcMaxCondElimination(i1,i2,i3, contPos,xn,tn,N,C,solBasis,solInhom );
s2=i2+i3*3;
condM=sparse(i3+i2,s2);
adderC=i2;

%%%%%%%%%%%%%%%%%%Kapazitätenbedingungen
for t=1:tn+1
    for i=1:N+1
        for j=1:N+1
            condM((t-1)*N1^2+(i-1)*N1+j,adderC+(t-1)*N1^2+(i-1)*N1+j)=1;
            if (i>1) && (j>1)%the input from source is independent of the heat flow
                condM((t-1)*N1^2+(i-1)*N1+j,i2+2*i3+(t-1)*N1^2+(i-1)*N1+j)=-C(i,j,1);    %Kapazitätenbedingung MIP
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
    V(1,i1+i)=1;
    for j=1:N+1
        V(1,adderC+(i)*N1+j)=1;
    end
end
A=[A;V];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Stateconstraints
A=[A; [solBasis sparse(i1,3*i3)]];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate b_L and b_U
[s1 s2]=size(TM);

b_U=zeros(i3+i2,1);
b_L=b_U;
b_L(1:i3)=-inf*ones(i3,1);
for t=1:tn+1
    for i=1:N+1
       b_U((t-1)*(N+1)^2+(i-1)*(N+1)+1)=C(i,1,1);
       b_U((t-1)*(N+1)^2+i)=C(1,i,1);
    end
end
b_U=[b_U;tempArcInhom];
b_L=[b_L;-inf*ones(s1+(tn+1)*(N)^2,1)];
b_U=[b_U;zeros(N^2,1)];
b_U=[b_U;(Tmax+Schwell)*ones((tn)*(N)^2,1)];
b_U=[b_U;0];
b_L=[b_L;0];
%state constraints for B
b_L=[b_L; -solInhom];
b_U=[b_U; inf*ones(i1,1)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%construct c
n=size(A,2);
c=zeros(n,1);
c(1:i2)=solBasis'*ones(tn1*(xn+1)^2,1);
c(n-i3+1:n)=-0.00001*ones(i3,1);

end