% function [ret,cut,sense,rhs] = cpxcb_USERCUT(x,f,xstat,Prob)
%
% CPLEX MIP User cut callback.
%
% The User Cut Callback is enabled by setting callback(15)=1 in the call to
% cplex.m, or Prob.MIP.callback(15)=1 if using tomRun('cplex',...)
%
% This callback is called by CPLEX during MIP branch & cut for every 
% node that has an LP optimal solution with objective value below the 
% cutoff and is integer infeasible. CPLEX also calls the callback when 
% comparing an integer feasible solution, including one provided by a 
% MIP start before any nodes exist, against lazy constraints.
%
% The callback routine can add globally valid cuts to the LP subproblem. A cut
% is a constraint of the following form:
%
%    c1*x(1) + c2*x(2) + ... + cn*x(n) <?>  rhs
% 
% where <?> is exactly one of the relations <=, >= or = and rhs is a scalar
% right hand side limit. 
%
% Calling syntax:
%
% function [ret,cut,sense,rhs] = cpxcb_USERCUT(x,f,xstat,Prob)
%  
% cpxcb_USERCUT is called by the solver with four arguments:
%
%  x     - The new integer solution
%  f     - The objective value at x
%  xstat - Feasibility information for each variable.  
%  Prob  - The TOMLAB problem structure
%
%  xstat elements may have values 0,1 or 2 with the following
%  interpretation for each variable in the current subproblem, 
%
%  0     - Variable is integer-valued
%  1     - Variable is not integer-valued
%  2     - Variable is "implied integer feasible"; it may have a fractional 
%          value in the current solution, but it will take on an integer 
%          value when all integer variables still in the problem have 
%          integer values. It should not be branched upon.
%
%
% cpxcb_USERCUT should return four values, 
%
% ret   - Should return one of the following scalar values:
%
%  0    - Continue optimization, use cuts as added (default)
%  1    - Terminate optimization
%         Any other return value will be interpreted as 0.
%  2    - Use cuts as added (same as 0)
%  3    - Exit cut loop and move to branching
%
% cut   - A sparse or dense matrix with cuts to add to the current LP
%         subproblem. This matrix must be a m*n sparse or dense double matrix,
%         where m is the number of cuts to add and n is the number of variables
%         (equal to length(x)). If no cuts are to be added, return [] (empty).
%
% sense - If the 'cut' return argument is nonempty, this array should be a m*1 
%         character array indicating the sense of each cut. Allowed values are 
%         L, <, G, >, E, =, for 'less than', 'greater than' and 'equal' respectively. 
%         If more than one cut is specified, simply stack the values, e.g.
%         'LEE'.
%
% rhs   - If the 'cut' return argument is nonempty, this array should be a m*1
%         array of right hand side values for each cut. Cuts are single-sided
%         constraints, thus only one rhs value exists for each cut.
%
% If modifying this file, it is recommended to make a copy of it which
% is placed before the original file in the MATLAB path.
%

% Anders Goran, Tomlab Optimization Inc., E-mail, tomlab@tomopt.com
% Copyright (c) 2002-2016 by Tomlab Optimization Inc., $Release, 12.6.1$
% Written Aug 1, 2007.  Last modified Jul 19, 2016.

function [ret,cut,sense,rhs] = cpxcb_USERCUT(x,f,xstat,wherefrom,Prob)
global solBasis;
global solInhom;
global xn;
global i2;
global tn;
global Aext;
wherefrom
% NOTE, Important that empty outputs are returned if no cuts are added. 
% if xstat==5
%     cut = [];  rhs = [];
%     state=solBasis*x(1:i2);
%     state=state+solInhom;
%     %Verbesserung:für jeden Zeitschritt i1 Bedingungen hinzufügen
%     for i=1:tn+1
%         [ASorted, AIdx] = sort(state((i-1)*(xn+1)^2+1:i*(xn+1)^2));
%          for k=1:1
%              addI=AIdx(k);
%              if state((i-1)*(xn+1)^2+addI)<-0.0001;
%                    cut=[cut; Aext((i-1)*(xn+1)^2+addI,:)];
%                    rhs=[rhs; -solInhom((i-1)*(xn+1)^2+addI)];
%              end
%          end
%     end
%     sense=repmat('>',length(rhs));?
% else
   cut=Aext;
   rhs=-solInhom;
   sense=repmat('>',length(rhs));
%end
% Insert code here to set values in cut, sense, rhs. 

% Return value, 0 to continue, 1 to terminate optimization. 
ret = 0;

% MODIFICATION LOG:
%
% 070801 ango Wrote function
% 131001 ango Updated to support latest CPLEX functionality


% Support routine: where are we called from? 
function w = cpxWhereFrom(wherefrom)
switch(wherefrom)
   case 101, w='MIP';
   case 102, w='MIP Branch';
   case 103, w='MIP Node';
   case 104, w='MIP Heuristic';
   case 105, w='MIP Solve';
   case 106, w='MIP Cut Loop';
   case 107, w='MIP Probe';
   case 108, w='MIP Fraccut';
   case 109, w='MIP Disjcut';
   case 110, w='MIP Flowmir';
   case 111, w='MIP Incumbent Nodesoln';
   case 112, w='MIP Deletenode';
   case 113, w='MIP Branch Nosoln';
   case 114, w='MIP Cut Last';
   case 115, w='MIP Cut Feas';
   case 116, w='MIP Cut Unbd';
   case 117, w='MIP Incumbent Heursoln';
   case 118, w='MIP Incumbent Usersoln';
   otherwise, w='Unknown';
end
   