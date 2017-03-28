%clear all;
%profile on;
%Name='optControl';
clear all;
itas=12:1:12;
count=1;
global xn;
global i2;
global tn;
%duras=zeros(length(itas),2);
for xn=itas
    tn=50;
    eliminate=1;

    video=1;
    saver=1;
    start=0;
    plotM2=0;
    MIP=1;
    CNFD=0;
    Tmax=10;
    Schwell=0.44; 
    %xn=12;
    time=3;
    dt=time/tn;
    dx=1/xn;
    dy=dx;
    a1=0.002;
    a2=0.000;
    d=0.01;
    N=25;
    i1=(xn+1)^2*(tn+1);
    i2=(tn)*N;
    i3=(tn+1)*(N+1)^2;
    G=zeros(N,2);
    r1=a1*dt/dx;
    r2=a2*dt/dx;
    R=d*dt/(dx*dx);
    cx=a1*dt/dx;
    cy=a2*dt/dx;
    Rx=a1*dx/d;
    Ry=a2*dx/d;
    %R=0;
    arctime=0.1;
    slowdown=dt/arctime;
    %slowdown=1;
    g=@(x,y) (x>=0.85)||(y>=0.85) || (x<=0.15) ||(y<=0.15) ;
    %g=@(x,y) 1;
    hL=10;
    hR=5;
    hU=10;
    hO=2;



    f=@(x,y,t) x>0.99 || y>0.99 || x<0.01 || y<0.01;
    g=@(x,y) 4*max(abs(x-0.5)^2,abs(y-0.5)^2);
    u=@(x,y,xc,yc) -0.1*exp(-10*(x-xc).^2-10*(y-yc).^2);
    for i=1:sqrt(N)
        for j=1:sqrt(N)
            G((i-1)*sqrt(N)+j,1)=1/(sqrt(N)+1)*j;
            G((i-1)*sqrt(N)+j,2)=1/(sqrt(N)+1)*i;
        end
    end

    C=constrC(0.4,0.15,N,5);

    duras(count,1)=runElimination(i1,i2,i3, R,xn,tn,dx,dt,N,u,G,a1,a2,r1,r2,C,g,f,Tmax,Schwell,slowdown,1,hL,hR,hU,hO,d,video);
    %duras(count,2)=runElimination(i1,i2,i3, R,xn,tn,dx,dt,N,u,G,a1,a2,r1,r2,C,g,f,Tmax,Schwell,slowdown,0,hL,hR,hU,hO,d,video);
    count=count+1;
end
%save('timedurastnvar');
%x_k=runElimination(i1,i2,i3, R,xn,tn,dx,dt,N,u,G,a1,a2,r1,r2,C,g,f,Tmax,Schwell,slowdown);
%profile viewer
%optVal=c(1:i1)'*x_k(1:i1)
%flow=x_k(i1+1:i1+i2);
%fctn=x_k(1:i1);


