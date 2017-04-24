%clear all;
%profile on;
%Name='optControl';


clear variables;
global xn;
global i2;
global tn;
global loader;
global usetime;
global time;
global contamination;
addpath('~/MIPDECO/Feuerprojekt/oFEM/projects/NEODYM/heat/');
addpath('~/MIPDECO/Feuerprojekt/oFEM/source/');
%%%%%controls
contamination=1
global plotM;
load('testTemp');

loader=1;
video=1;
eliminate=1;
saver=1;
scenario=5;
xn=25;
tn=40;
%%%%%setup
if loader
    load(sprintf('statexn%dtn%d.mat',xn,tn));
    load(sprintf('matrixData%d_%d.mat',xn,tn));
    xn1=xn+1;
    tn1=tn+1;
    if video
        dx=1/xn;
        dt=1/tn;
        x_k=x_k';
        %x_k=Result.x_k;
        x_kEli=x_k;
        fac=1;
        slowdown=dt/arctime;
        step=ceil(round(1/slowdown,2));
        if contamination
            plotVideoContamination
        else
            plotVideoVertex
        end
    end
else
    switch scenario
        case 1
            [C,G,N]=constrGraph(0.2,0.15,7);
            params=struct('xn',40,'tn',30,'Tmax',5,'Schwell',0.8,'slowdown',0.1,'N',N,'time',1,'u',@(x,xc) ofem.matrixarray(-30*exp(-50*dot(x-xc,x-xc,1))));
            [paramsControlled,paramsInhom]=PDEparams(1);
            arctime=params.time/30;
            usetime=params.time/20;
            
        case 2
            [C,G,N]=constrGraph(0.2,0.15,8);
            params=struct('xn',45:5:45,'tn',30:10:30,'Tmax',5,'Schwell',0.8,'slowdown',0.1,'N',N,'time',1,'u',@(x,xc) ofem.matrixarray(-25*exp(-30*dot(x-xc,x-xc,1))));
            [paramsControlled,paramsInhom]=PDEparams(2);
            arctime=params.time/30;
            usetime=params.time/20;
            
        case 3
            [C,G,N]=constrGraph(0.3,0.15,9);
            params=struct('xn',45,'tn',30,'Tmax',2,'Schwell',0.8,'slowdown',0.1,'N',N,'time',1,'u',@(x,xc) ofem.matrixarray(-30*exp(-50*dot(x-xc,x-xc,1))));
            [paramsControlled,paramsInhom]=PDEparams(3);
            arctime=params.time/30;
            usetime=params.time/15;
            
        case 4
            [C,G,N]=constrGraph(0.3,0.15,10);
            params=struct('xn',80,'tn',30,'Tmax',2,'Schwell',0.8,'slowdown',0.1,'N',N,'time',1,'u',@(x,xc) ofem.matrixarray(-10*exp(-50*dot(x-xc,x-xc,1))));
            [paramsControlled,paramsInhom]=PDEparams(4);
            arctime=params.time/30;
            usetime=params.time/15;
        case 5
            [C,G,N]=constrGraph(0.2,0.15,12);
            params=struct('xn',25:5:25,'tn',40:10:40,'Tmax',5,'Schwell',0.8,'slowdown',0.1,'N',N,'time',1,'u',@(x,xc) ofem.matrixarray(-25*exp(-30*dot(x-xc,x-xc,1))));
            [paramsControlled,paramsInhom]=PDEparams(5);
            arctime=params.time/30;
            usetime=params.time/20;
            
    end
    %___________________

    %contPos=[0.5,0.5];
    %T_xt=heatRobin(gridPoints,0:0.025:1,@(x) ofem.matrixarray(-15000*exp(-10*dot(x-contPos',x-contPos',1))),@(x) 100+0*x(1,1,:),struct('c_P',0.17,'rho',750,'k',2400,'h',0.1,'v',[0;0],'T_a',1000),struct('LL',[0,0],'UR',[1,1],'h_max',1/20));
    itas=22:2:22;
    count=1;

    time=params.time; %the time the water has an effect.
    for xn=params.xn
        for tn=params.tn
            save(sprintf('matrixData%d_%d.mat',xn,tn),'params','paramsControlled','paramsInhom','arctime','usetime','C','G','N');
            dt=time/tn;
            dx=1/xn;
            dy=dx;
            i1=(xn+1)^2*(tn+1);
            i2=(tn)*N;
            i3=(tn+1)*(N+1)^2;
            slowdown=dt/arctime;
            %
            g=@(x,y) 0.3+1.2*((x-1)^2 +(y-1)^2<0.4);
            if contamination
                u= @(x,y,xc,yc) -0.01*exp(-30*((x-xc)^2+(y-yc)^2));
            else
                u= @(x,y,xc,yc) -25*exp(-30*((x-xc)^2+(y-yc)^2));
            end
            
            %runVertex(i1,i2,i3,xn,tn,dx,dt,N,u,G,-0.01,-0.01,  C,g,params.Tmax,params.Schwell,0,0,0,0,0.02,slowdown);
            runEliminationVertex(i1,i2,i3,xn,tn,dx,dt,params.N,params.u,G,C,params.Tmax,params.Schwell,slowdown,1,video,paramsControlled,paramsInhom);
        end
    end
end
