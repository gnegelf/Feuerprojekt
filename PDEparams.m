function [ paramsControlled,paramsInhom,otherParams ] = PDEparams( scenario )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
switch scenario
    case 1
        paramsControlled=struct('c_P',1,'rho',1,'k',0.02,'h',2,'v',[-0.1;-0.1],'T_a',[0,0,0,0],'initialT',@(x) ofem.matrixarray(0.5+0*x(1,1,:)));
        paramsInhom=paramsControlled;
        paramsInhom.T_a=[0.5,0.5,5,5];
    case 2
        %paramsControlled=struct('c_P',1,'rho',1,'k',0.02,'h',0,'v',[-0.1;-0.1],'T_a',[0,0,0,0],'initialT',@(x) ofem.matrixarray(0.3+2*(x(1,1,:)>[0.5,0.5])));
        paramsControlled=struct('c_P',1,'rho',1,'k',0.02,'h',0,'v',[-0.01;-0.01],'T_a',[0,0,0,0],'initialT',@(x) ofem.matrixarray(0.3+1.2*((dot(x-[1;1],x-[1;1],1))<0.4)));
        paramsInhom=paramsControlled;
        paramsInhom.T_a=[0.4,0.4,0.4,0.4];
    case 3
        %paramsControlled=struct('c_P',1,'rho',1,'k',0.02,'h',0,'v',[-0.1;-0.1],'T_a',[0,0,0,0],'initialT',@(x) ofem.matrixarray(0.3+2*(x(1,1,:)>[0.5,0.5])));
        paramsControlled=struct('c_P',1,'rho',1,'k',0.02,'h',1,'v',[-0.08;-0.0],'T_a',[0,0,0,0],'initialT',@(x) ofem.matrixarray(0.3+1.2*(x(1,1,:)>0.7)));
        paramsInhom=paramsControlled;
        paramsInhom.T_a=[0.3,0.3,2.0,0.3];
    case 4
        paramsControlled=struct('c_P',1,'rho',1,'k',0.02,'h',1,'v',[-0.03;-0.0],'T_a',[0,0,0,0],'initialT',@(x) ofem.matrixarray(0.3+1.2*(x(1,1,:)>0.3 & x(1,1,:)<0.7 )));
        paramsInhom=paramsControlled;
        paramsInhom.T_a=[0.3,0.3,0.3,0.3];
    case 11
        paramsControlled=struct('c_P',1,'rho',1,'k',0.007,'h',[1,0,1,0],'v',[-0.1;0.0],'T_a',[0,0,0,0],'initialT',@(x) ofem.matrixarray(0.5+0*x(1,1,:)));
        paramsInhom=paramsControlled;
        paramsInhom.T_a=[0,0,5,0];
        
end

end

