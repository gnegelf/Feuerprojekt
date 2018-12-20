function [ ] = runEliminationVertex(i1,i2,i3,xn,tn,dx,dt,N,u,G,C,Tmax,Schwell,slowdown,eliminate,video, p1,p2)
%%%This File runs the fire Model with linear heat equation, without
%%%eliminating the PDE it is a cleaned up version of the previous trials
%Lineare PDE, lineare Abh�ngigkeit vom Controller, Zus�tzliche Netzwerk
%Bedingungen
tic
fac=1;
C(:,:,1)=C(:,:,1)*fac;
step=ceil(round(1/slowdown,2));
global contamination;
%[A]= constrACNFD( i1,i2,i3,R,xn,tn,dx,N,u,G,a1,a2,r1,r2,TM,ceil(1/slowdown),C,Tmax,r1,r2,a1*dx/d,a2*dx/d);
global Aext;
global solInhom;
global solBasis;
global loader;
global usetime;
global scenario;
%C(:,:,1)=C(:,:,1)/step;
if loader
    load(sprintf('data/statexn%dtn%d_%d.mat',xn,tn,scenario));
    load(sprintf('data/matrixData%d_%d_%d.mat',xn,tn,scenario));
else
    if contamination
        usetime=dt;
         [A,b_U,b_L,c,solInhom,solBasis,Aext,b_Uext,b_Lext,A_ext2,b_Uext2,b_Lext2]= constrContamination(i1,i2,i3 ,xn,tn,dx,dt,N,u,G,step,C,Tmax,Schwell,p1,p2);
         [m,n]=size(A);
         intVarN=N;
    else        
        [A,b_U,b_L,c,solInhom,solBasis,Aext,b_Uext,b_Lext,A_ext2,b_Uext2,b_Lext2]= constrMIPEliminationRobinVertexFECONT(i1,i2,i3 ,xn,tn,dx,dt,N,u,G,step,C,Tmax,Schwell,p1,p2);
        [m,n]=size(A);
        intVarN=i2;
    end
    %C(:,:,1)=C(:,:,1)*step;


    
    %IntVars=[];
    %Prob      = lpAssign(c, A, b_L, b_U, x_L, x_U, x_0, Name, setupFile, ...
    %                      nProblem, IntVars, VarWeight, KNAPSACK, fIP, xIP, ...
    %                      f_Low, x_min, x_max, f_opt, x_opt);


    contVarN=n-intVarN;
    Afull=[A;Aext];
    %Afull=full(Afull);
    %A=full(A);
    %Aext=full(Aext);

    b_Ufull=[b_U;b_Uext];
    b_Lfull=[b_L;b_Lext];

    A2=[A;A_ext2];
    b_U2=[b_U;b_Uext2];
    b_L2=[b_L;b_Lext2];
    if ~contamination
        saveName=sprintf('~/MIPDECO/Feuerprojekt/data/matrixData%d_%d_%d.mat',xn,tn,scenario);
        save(saveName,'A','b_U','b_L','c','solInhom','solBasis','-append');
        saveName=sprintf('~/python/data/feuerData%d_%d_%d.mat',xn,tn,scenario);
        save(saveName,'A2','-v7.3','Aext','b_Uext','b_Lext','xn','tn','intVarN','contVarN','b_U2','b_L2','c','Afull','b_Ufull','b_Lfull')
    else
        saveName=sprintf('~/MIPDECO/Feuerprojekt/data/contaMatrixData%d_%d_%d.mat',xn,tn,scenario);
        save(saveName,'A','b_U','b_L','c','solInhom','solBasis','-append');
        saveName=sprintf('~/python/data/contaData%d_%d_%d.mat',xn,tn,scenario);
        save(saveName,'A2','-v7.3','Aext','b_Uext','b_Lext','xn','tn','intVarN','contVarN','b_U2','b_L2','c','Afull','b_Ufull','b_Lfull')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%count=100000;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Result = tomRun('cplex', Prob2, 4);
if loader
    x_k=x_k22;
    x_k=x_k';
    %x_k=Result.x_k;
    x_kEli=x_k;
    plotVideoVertex
end
end