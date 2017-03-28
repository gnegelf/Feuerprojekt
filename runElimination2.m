function [ duration] = runElimination(i1,i2,i3, R,xn,tn,dx,dt,N,u,G,a1,a2,r1,r2,C,g,f,Tmax,Schwell,slowdown,eliminate,hL,hR,hU,hO,d,video)
%%%This File runs the fire Model with linear heat equation, without
%%%eliminating the PDE it is a cleaned up version of the previous trials
%Lineare PDE, lineare Abh�ngigkeit vom Controller, Zus�tzliche Netzwerk
%Bedingungen
tic
fac=1;
C(:,:,1)=C(:,:,1)*fac;
step=ceil(round(1/slowdown,2));

%[A]= constrACNFD( i1,i2,i3,R,xn,tn,dx,N,u,G,a1,a2,r1,r2,TM,ceil(1/slowdown),C,Tmax,r1,r2,a1*dx/d,a2*dx/d);
global Aext;
global solInhom;
global solBasis;
%C(:,:,1)=C(:,:,1)/step;
if eliminate 
    [A,b_U,b_L,c,solInhom,solBasis,Aext,b_Uext,b_Lext]= constrMIPEliminationRobin( i1,i2,i3,R,xn,tn,dx,dt,N,u,G,a1,a2,step,C,g,f,Tmax,Schwell,hL,hR,hU,hO,d);
    %C(:,:,1)=C(:,:,1)*step;

    
    [m,n]=size(A);
    xn1=xn+1;
    tn1=tn+1;
    x_L=zeros(n,1);
    x_U=100000*ones(n,1);
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
    intVarN=length(IntVars);
    contVarN=n-intVarN;
    Afull=[A;Aext];
    %Afull=full(Afull);
    %A=full(A);
    %Aext=full(Aext);
    b_Ufull=[b_U;b_Uext];
    b_Lfull=[b_L;b_Lext];
    save('feuerData.mat','A','Aext','b_Uext','b_Lext','xn','tn','intVarN','contVarN','b_U','b_L','c','Afull','b_Ufull','b_Lfull')
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
    %while (count<N*5)
        count=count+1;
        [m,n]=size(A);
        Prob      = mipAssign(c, A, b_L, b_U, x_L, x_U, x_0, 'contamination', setupFile, ...
                              nProblem, IntVars, VarWeight, KNAPSACK, fIP, xIP, ...
                              f_Low, x_min, x_max, f_opt, x_opt);
        %Prob2= preSolve(Prob);
        %Prob.optParam.IterPrint = 0; % Set to 1 to see iterations.
        %Prob.cpxControl.LPMETHOD = 1;
        Prob.Solver.Alg = 5;   % Depth First, then Breadth search
        Prob.PriLev =3;
        Prob.optParam.IterPrint = 3;
        Prob.LGO.tlimit=5000;
        Prob.MIP.callback = zeros(16,1);
        Prob.MIP.callback(16) = 1;
        % Calling driver routine tomRun to run the solver.
        % The 1 sets the print level after optimization.
        Result = tomRun('cplex', Prob, 4);
%         x_k=Result.x_k;
%         state=solBasis*x_k(1:i2);
%         state=state+solInhom;
%         Aadd=[];
%         b_Uadd=[];
%         b_Ladd=[];
%         %Verbesserung:f�r jeden Zeitschritt i1 Bedingungen hinzuf�gen
%         for i=1:tn+1
%             [ASorted, AIdx] = sort(state((i-1)*(xn+1)^2+1:i*(xn+1)^2));
%              for k=1:1
%                  addI=AIdx(k);
%                  if state((i-1)*(xn+1)^2+addI)<-0.0001;
%                        Aadd=[Aadd; Aext((i-1)*(xn+1)^2+addI,:)];
%                        b_Ladd=[b_Ladd; -solInhom((i-1)*(xn+1)^2+addI)];
%                        b_Uadd=[b_Uadd; inf];
%                  end
%              end
%         end
%         A=[A; Aadd];
%         b_U=[b_U; b_Uadd];
%         b_L=[b_L;b_Ladd];
%         if sum(state<-0.0001)==0
%             break;
%         end
        count
    %end

    duration=toc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Result = tomRun('cplex', Prob2, 4);
    x_k=Result.x_k;
    x_kEli=x_k;
    plotVideo
else
    [A,b_U,b_L,c]=constrANoElimination(i1,i2,i3, R,xn,tn,dx,dt,N,u,G,a1,a2,step,C,g,f,Tmax,Schwell,hL,hR,hU,hO,d);
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

    Prob      = mipAssign(c, A, b_L, b_U, x_L, x_U, x_0, 'contamination', setupFile, ...
                          nProblem, IntVars, VarWeight, KNAPSACK, fIP, xIP, ...
                          f_Low, x_min, x_max, f_opt, x_opt);
    

    %Prob2= preSolve(Prob);
    %Prob.optParam.IterPrint = 0; % Set to 1 to see iterations.
    %Prob.cpxControl.LPMETHOD = 1;
    Prob.Solver.Alg = 5;   % Depth First, then Breadth search
    Prob.PriLev =3;
    Prob.optParam.IterPrint = 3;
    Prob.LGO.tlimit=5000;
    % Calling driver routine tomRun to run the solver.
    % The 1 sets the print level after optimization.
    Result = tomRun('cplex', Prob, 4);
    x_k=Result.x_k((xn+3)^2*tn1+1:length(Result.x_k));
    
    duration=toc
    x_kNoEli=x_k;
    plotVideo;
end
end