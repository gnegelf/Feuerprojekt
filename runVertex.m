function [ ] = runVertex(i1,i2,i3,xn,tn,dx,dt,N,u,G,a1,a2,C,g,Tmax,Schwell,hL,hR,hU,hO,D,slowdown,loadedSolution)
%%%This File runs the fire Model with linear heat equation, without
%%%eliminating the PDE it is a cleaned up version of the previous trials
%Lineare PDE, lineare Abh�ngigkeit vom Controller, Zus�tzliche Netzwerk
%Bedingungen
tic
fac=1;
C(:,:,1)=C(:,:,1)*fac;
step=ceil(round(1/slowdown,2));
global usetime;
global scenario;

if ~loadedSolution
    [A,b_U,b_L,c]=  constrANoElimination(i1,i2,i3, xn,tn,dx,dt,N,u,G,a1,a2,step,C,g,Tmax,Schwell,hL,hR,hU,hO,D);
    [m,n]=size(A);
    intVarN=i2;
    contVarN=n-intVarN;
    saveName=sprintf('matrixData%d_%d_%d.mat',xn,tn,scenario);
    save(saveName,'A','b_U','b_L','c','-append');
    saveName=sprintf('~/python/feuerDataNoElimination%d_%d.mat',xn,tn);
    save(saveName,'A','b_U','b_L','xn','tn','intVarN','contVarN','c');
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%count=100000;
video=1;
loadName=sprintf('stateNoElixn%dtn%d.mat',xn,tn);
load(loadName);
xn1=xn+1;
tn1=tn+1;
dx=1/xn;
dt=1/tn;
x_k=x_k';
%x_k=Result.x_k;
x_kEli=x_k;
fac=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotVideoVertexNoEli;
end
end