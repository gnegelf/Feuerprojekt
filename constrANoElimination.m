function [ A,b_U,b_L,c ] = constrANoElimination(i1,i2,i3, xn,tn,dx,dt,N,u,contPos,a1,a2,step,C,g,Tmax,Schwell,hL,hR,hU,hO,D)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
i4=(xn+3)^2*(tn+1);
[TM]=calcMaxCond(i4,i2,i3, contPos,xn,tn,N,Tmax );
xn1=xn+3;
tn1=tn+1;
N1=N+1;
%A=speye(i4,i4+i2+i3*3);
R=D*dt/dx^2;
%discretization of Laplace
A=sparse(xn1^2*(tn+1),xn1^2*(tn+1)+2*i2+i3);
A(1:xn1^2,1:xn1^2)=speye(xn1^2,xn1^2);
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
A(1:xn1^2,1:xn1^2)=speye(xn1^2,xn1^2);
Aadd=kron(speye(tn),difM(:,1:xn1^2));
Aadd2=kron(speye(tn),difM(:,xn1^2+1:xn1^2*2));
A(xn1^2+1:xn1^2*(tn+1),1:xn1^2*(tn))=Aadd;
A(xn1^2+1:xn1^2*(tn+1),xn1^2+1:xn1^2*(tn+1))=A(xn1^2+1:xn1^2*(tn+1),xn1^2+1:xn1^2*(tn+1))+Aadd2;
global usetime;
%controlM=[controlM;zeros(N,N)];
%difMext=[difM controlM];

%Kontrolleinfluss
for t=2:tn+1
    %A((t-1)*(xn1^2)+1:t*(xn1^2),(t-2)*(xn1^2)+1:t*xn1^2)=difM;
    for i=2:xn+2
        for j=2:xn+2
            for k=1:N
                for tt=0:min(usetime,tn+1-t)
                    if abs(u((i-2)*dx,(j-2)*dx,contPos(k,1),contPos(k,2)))>0.0001
                        A((t+tt-1)*xn1^2+(i-1)*xn1+j,tn1*xn1^2+(t-2)*N+k)=-u((i-2)*dx,(j-2)*dx,contPos(k,1),contPos(k,2));%vorher: u((i-1)*dx,(j-1)*dx,contPos(k,1),contPos(k,2));
                    end
                end
           end
        end
    end
end
global plotM;
b=zeros(i4,1);
for t=1:tn+1
    if t==1
        for k=2:xn+2
            for j=2:xn+2

                b((t-1)*(xn1)^2+(j-1)*(xn1)+k)=plotM(j-1,k-1,1);
                %b((t-1)*(xn1)^2+(j-1)*(xn1)+k)=g(dx*(j-2),dx*(k-2));
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


[s1,s2]=size(A); %#ok<ASGLU>
condM=sparse(i3+2*i2,s2);
adderC=i4+i2;

%Kapazit�tenbedingung MIP
for t=1:tn+1
    for i=1:N+1
        for j=1:N+1
            condM((t-1)*N1^2+(i-1)*N1+j,adderC+(t-1)*N1^2+(i-1)*N1+j)=1;
            if (i>1) && (j>1)
                %condM((t-1)*N1^2+(i-1)*N1+j,i4+i2+2*i3+(t-1)*N1^2+(i-1)*N1+j)=-C(i,j,1);
                for cc=0:C(i,j,2)-1
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


for t=2:tn+1
    for i=1:N
        condM(i3+(t-2)*N+i,i4+(t-2)*N+i)=1; %set coefficient of variable corresponding to wateroutput at controller i at time t-1
        for j=1:N+1
            condM(i3+(t-2)*N+i,adderC+(t-1)*N1^2+(i)*N1+j)=1;
            if t>step*C(i+1,j,2)
                condM(i3+(t-2)*N+i,adderC+(t-step*C(i+1,j,2)-1)*N1^2+(j-1)*N1+i+1)=-1;
            end
        end
    end   
end


for t=1:tn
    for i=1:N
        condM(i3+i2+(t-1)*N+i,i4+(t-1)*N+i)=1;
        condM(i3+i2+(t-1)*N+i,i4+i3+i2+(t-1)*N+i)=-1;
    end
end
count=0;
for t=2:tn+1
    for i=1:N
        for j=1:N
            if C(j,i,1) > 0
                count=count+1;
                condM(i3+2*i2+count,i4+i2+(t-1)*(N+1)^2+j*(N+1)+i+1)=1;
                condM(i3+2*i2+count,i4+i3+i2+(t-2)*N+i)=-1;
            end
            
        end
    end
end

A=[A;condM];
A=[A;TM];
%A=[A;zeros(size(TM))];
V=zeros(1,s2);
%start flow = zero condition
for i=1:N
    V(1,i4+i)=1;
    for j=1:N+1
        V(1,adderC+(i)*N1+j)=1;
    end
end
A=[A;V];
[s1 s2]=size(TM);


b_U=[b;zeros(i3+i2+i2+count,1)];%LGS für Temperatur i4 Bed
b_L=b_U;

%Kapazitätenbedingungen i3 Bed.
for t=1:tn+1
    for i=1:N+1
        for j=1:N+1
           b_U(i4+(t-1)*(N+1)^2+(i-1)*(N+1)+j)=C(i,j,1);
        end
    end
end
b_L(i4+i3+i2+1:i4+i3+i2+i2+count)=-1*ones(i2+count,1);
b_U=[b_U;zeros(s1,1)+Tmax+Schwell];
b_L=[b_L;-10000*ones(s1,1)];

b_U=[b_U;0];
b_L=[b_L;0];
%

[m,n]=size(A);
c=zeros(n,1);
c(1:i4)=ones(i4,1);
for t=1:tn+1
    for x=1:xn1
        for y=1:xn1
            if x==1 || y==1 || x==xn1 || y==xn1
                c((t-1)*xn1^2+(x-1)*xn1+y)=0; 
            end
        end
    end
end
c(n-i3+1:n)=-0.00001*ones(i3,1);
end