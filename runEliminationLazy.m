function [ duration] = runEliminationLazy(i1,i2,i3, R,xn,tn,dx,dt,N,u,G,a1,a2,r1,r2,C,g,f,Tmax,Schwell,slowdown,eliminate,hL,hR,hU,hO,d,video)
%%%This File runs the fire Model with linear heat equation, without
%%%eliminating the PDE it is a cleaned up version of the previous trials
%Lineare PDE, lineare Abhängigkeit vom Controller, Zusätzliche Netzwerk
%Bedingungen
tic
fac=1;
C(:,:,1)=C(:,:,1)*fac;
step=ceil(round(1/slowdown,2));

%[A]= constrACNFD( i1,i2,i3,R,xn,tn,dx,N,u,G,a1,a2,r1,r2,TM,ceil(1/slowdown),C,Tmax,r1,r2,a1*dx/d,a2*dx/d);

%C(:,:,1)=C(:,:,1)/step;
[A,b_U,b_L,c,solInhom,solBasis,Aext]= constrMIPEliminationRobin( i1,i2,i3,R,xn,tn,dx,dt,N,u,G,a1,a2,step,C,g,f,Tmax,Schwell,hL,hR,hU,hO,d);
%C(:,:,1)=C(:,:,1)*step;


[m,n]=size(A);
xn1=xn+1;
tn1=tn+1;
x_L=zeros(n,1);
x_U=inf*ones(n,1);
x_U(n-i3+1:n)=ones(i3,1);
IntVars=n-i3+1:n;
%IntVars=[];
x_min   = x_L; 
x_max  = x_U; 
x_0=[];
%Prob      = lpAssign(c, A, b_L, b_U, x_L, x_U, x_0, Name, setupFile, ...
%                      nProblem, IntVars, VarWeight, KNAPSACK, fIP, xIP, ...
%                      f_Low, x_min, x_max, f_opt, x_opt);
f_Low  = -1E7; 
f_opt   = [];
nProblem  = []; % Problem number not used
fIP       = []; % Do not use any prior knowledge
xIP       = []; % Do not use any prior knowledge
setupFile = []; % Just define the Prob structure, not any permanent setup file
x_opt     = []; % The optimal integer solution is not known
VarWeight = []; % No variable priorities, largest fractional part will be used
KNAPSACK  = []; 

count=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%count=100000;

% Prob      = mipAssign(c, A, b_L, b_U, x_L, x_U, x_0, 'contamination', setupFile, ...
%                           nProblem, IntVars, VarWeight, KNAPSACK, fIP, xIP, ...
%                           f_Low, x_min, x_max, f_opt, x_opt);
% Prob2=Prob;
%     %Prob.optParam.IterPrint = 0; % Set to 1 to see iterations.
%     %Prob.cpxControl.LPMETHOD = 1;
%     Prob2.Solver.Alg = 5;   % Depth First, then Breadth search
%     Prob2.PriLev =3;
%     Prob2.optParam.IterPrint = 3;
%     % Calling driver routine tomRun to run the solver.
%     % The 1 sets the print level after optimization.
%     Result = tomRun('cplex', Prob2, 4);                      
Prob      = mipAssign(c, A, b_L, b_U, x_L, x_U, x_0, 'contamination', setupFile, ...
                      nProblem, IntVars, VarWeight, KNAPSACK, fIP, xIP, ...
                      f_Low, x_min, x_max, f_opt, x_opt);

%Prob2= preSolve(Prob);

%Prob.optParam.IterPrint = 0; % Set to 1 to see iterations.
%Prob.cpxControl.LPMETHOD = 1;
Prob.Solver.Alg = 5;   % Depth First, then Breadth search
Prob.PriLev =3;
Prob.optParam.IterPrint = 3;
%Prob2.LGO.tlimit=5000;
lazy = struct( 'A', Aext', 'rhs', -solInhom, 'sense', repmat('>',length(solInhom),1));
Prob.CPLEX.lazycon = lazy;

% Calling driver routine tomRun to run the solver.
% The 1 sets the print level after optimization.
Result = tomRun('cplex', Prob, 4);
x_k=Result.x_k;

duration=toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Result = tomRun('cplex', Prob2, 4);
x_k=Result.x_k;
x_kEli=x_k;
plotVideo

end