function [ TMr,b ] = calcMaxCondEliminationVertex (i1,i2,i3, contPos,xn,tn,N,C,solBasis,solInhom,Tmax)

%This function sets the matrix entries of A s.t. a control at a vertex is
%not possible if the temparature exceeds Tschwell

contPos=contPos*(xn);
s2=i1+2*i2+i3;
TM=sparse(N*(tn+1),s2);
%b1=zeros(1,N+1);

for t=1:tn
    for i=2:N+1
        x=contPos(i-1,1);
        y=contPos(i-1,2);
        xL=floor(x);
        xU=ceil(x);
        yL=floor(y);
        yU=ceil(y);
        if abs((sqrt((xL-x)^2+(yL-y)^2)+sqrt((xL-x)^2+(yU-y)^2)+sqrt((xU-x)^2+(yL-y)^2)+sqrt((xU-x)^2+(yU-y)^2)))<0.0001
            TM((t-1)*N+i-1,(t-1)*(xn+1)^2+xL*(xn+1)+yL+1)=1;
        else
            factor=1/(sqrt((xL-x)^2+(yL-y)^2)+sqrt((xL-x)^2+(yU-y)^2)+sqrt((xU-x)^2+(yL-y)^2)+sqrt((xU-x)^2+(yU-y)^2));
            TM((t-1)*N+i-1,(t-1)*(xn+1)^2+xL*(xn+1)+yL+1)= TM((t-1)*N+i-1,(t-1)*(xn+1)^2+xL*(xn+1)+yL+1)+factor*sqrt((xL-x)^2+(yL-y)^2);%weight such that they sum up to one
            TM((t-1)*N+i-1,(t-1)*(xn+1)^2+xL*(xn+1)+yU+1)= TM((t-1)*N+i-1,(t-1)*(xn+1)^2+xL*(xn+1)+yU+1)+factor*sqrt((xL-x)^2+(yU-y)^2);
            TM((t-1)*N+i-1,(t-1)*(xn+1)^2+xU*(xn+1)+yL+1)= TM((t-1)*N+i-1,(t-1)*(xn+1)^2+xU*(xn+1)+yL+1)+factor*sqrt((xU-x)^2+(yL-y)^2);
            TM((t-1)*N+i-1,(t-1)*(xn+1)^2+xU*(xn+1)+yU+1)= TM((t-1)*N+i-1,(t-1)*(xn+1)^2+xU*(xn+1)+yU+1)+factor*sqrt((xU-x)^2+(yU-y)^2);
        end
        %TM((t-1)*N+i-1,i1+i2+i3+(t-1)*N+i-1)=Tmax; wrong matrix!!!
    end
end
% 
% for t=1:tn
%     for i=2:N+1
%         x=contPos(i-1,1);
%         y=contPos(i-1,2);
%         %xL=floor(x);
%         %xU=ceil(x);
%         %yL=floor(y);
%         %yU=ceil(y);
%         yL=floor(x);
%         yU=ceil(x);
%         xL=floor(y);
%         xU=ceil(y);
%         if abs((sqrt((xL-y)^2+(yL-x)^2)+sqrt((xL-y)^2+(yU-x)^2)+sqrt((xU-y)^2+(yL-x)^2)+sqrt((xU-y)^2+(yU-x)^2)))<0.0001
%             TM((t-1)*N+i-1,(t)*(xn+1)^2+xL*(xn+1)+yL+1)=1;
%         else
%             factor=1/(sqrt((xL-y)^2+(yL-x)^2)+sqrt((xL-y)^2+(yU-x)^2)+sqrt((xU-y)^2+(yL-x)^2)+sqrt((xU-y)^2+(yU-x)^2));
%             TM((t-1)*N+i-1,(t)*(xn+1)^2+xL*(xn+1)+yL+1)= TM((t-1)*N+i-1,(t)*(xn+1)^2+xL*(xn+1)+yL+1)+factor*sqrt((xL-y)^2+(yL-x)^2);%weight such that they sum up to one
%             TM((t-1)*N+i-1,(t)*(xn+1)^2+xL*(xn+1)+yU+1)= TM((t-1)*N+i-1,(t)*(xn+1)^2+xL*(xn+1)+yU+1)+factor*sqrt((xL-y)^2+(yU-x)^2);
%             TM((t-1)*N+i-1,(t)*(xn+1)^2+xU*(xn+1)+yL+1)= TM((t-1)*N+i-1,(t)*(xn+1)^2+xU*(xn+1)+yL+1)+factor*sqrt((xU-y)^2+(yL-x)^2);
%             TM((t-1)*N+i-1,(t)*(xn+1)^2+xU*(xn+1)+yU+1)= TM((t-1)*N+i-1,(t)*(xn+1)^2+xU*(xn+1)+yU+1)+factor*sqrt((xU-y)^2+(yU-x)^2);
%         end
%         %TM((t-1)*N+i-1,i1+i2+i3+(t-1)*N+i-1)=Tmax; wrong matrix!!!
%     end
% end

solExt=sparse(s2,2*i2+i3);
solExt(1:i1,1:i2)=solBasis;
inhomExt=zeros(s2,1);
inhomExt(1:i1)=solInhom;
TMr=TM*solExt;
for t=1:tn
    for i=2:N+1
        TMr((t-1)*N+i-1,i2+i3+(t-1)*N+i-1)=Tmax;
    end
end
b=TM*inhomExt;
end
