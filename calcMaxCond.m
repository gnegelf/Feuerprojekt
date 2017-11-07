function [ TM,b ] = calcMaxCond (i1,i2,i3, contPos,xn,tn,N,Tmax )
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here
s2=i1+2*i2+i3;
contPos=contPos*(xn);
TM=sparse(N*(tn+1),s2);
%b1=zeros(1,N+1);
b=zeros(N*(tn+1),s2);
for t=1:tn
    for i=2:N+1
        x=contPos(i-1,1);
        y=contPos(i-1,2);
        xL=floor(x);
        xU=ceil(x);
        yL=floor(y);
        yU=ceil(y);
        TM((t-1)*N+i-1,i1+i2+i3+(t-1)*N+i-1)=Tmax;
        if abs((sqrt((xL-x)^2+(yL-y)^2)+sqrt((xL-x)^2+(yU-y)^2)+sqrt((xU-x)^2+(yL-y)^2)+sqrt((xU-x)^2+(yU-y)^2)))<0.0001
            TM((t-1)*N+i-1,(t-1)*(xn+3)^2+xL*(xn+3)+yL+1)=1;
        else
            factor=1/(sqrt((xL-x)^2+(yL-y)^2)+sqrt((xL-x)^2+(yU-y)^2)+sqrt((xU-x)^2+(yL-y)^2)+sqrt((xU-x)^2+(yU-y)^2));
            TM((t-1)*N+i-1,(t-1)*(xn+3)^2+xL*(xn+3)+yL+1)= TM((t-1)*N+i-1,(t-1)*(xn+3)^2+xL*(xn+3)+yL+1)+factor*sqrt((xL-x)^2+(yL-y)^2);%weight such that they sum up to one
            TM((t-1)*N+i-1,(t-1)*(xn+3)^2+xL*(xn+3)+yU+1)= TM((t-1)*N+i-1,(t-1)*(xn+3)^2+xL*(xn+3)+yU+1)+factor*sqrt((xL-x)^2+(yU-y)^2);
            TM((t-1)*N+i-1,(t-1)*(xn+3)^2+xU*(xn+3)+yL+1)= TM((t-1)*N+i-1,(t-1)*(xn+3)^2+xU*(xn+3)+yL+1)+factor*sqrt((xU-x)^2+(yL-y)^2);
            TM((t-1)*N+i-1,(t-1)*(xn+3)^2+xU*(xn+3)+yU+1)= TM((t-1)*N+i-1,(t-1)*(xn+3)^2+xU*(xn+3)+yU+1)+factor*sqrt((xU-x)^2+(yU-y)^2);
        end
        %TM((t-1)*N+i-1,i1+i2+i3+(t-1)*N+i-1)=Tmax; wrong matrix!!!
    end
end

% 
% oldCode
% s2=ii1+ii2+3*ii3;
% xn1=xn+3;
% contPos=contPos*(xn);
% b2=zeros(1,(N+1)^2);
% b1=zeros(1,(xn1)^2);
% T1=[];
% T2=[];
% for c1=2:N+1
%     for c2=2:N+1
%         if C(c1,c2,1)>0
%             x1=contPos(c1-1,1);
%             x2=contPos(c2-1,1);
%             y1=contPos(c1-1,2);
%             y2=contPos(c2-1,2);
%             dx=x2-x1;
%             dy=y2-y1;
%             if not(round(dx,5)==0)
%                m=dy/dx; 
%                 for i=round(x1):sign(x2-x1):round(x2)
% 
%                     j1=ceil(m*(i-x1)+y1);
%                     j2=floor(m*(i-x1)+y1);
%                     b1(1+(j1+1)*(xn1)+i+1)=1-(j1-(m*(i-x1)+y1));
%                     b1(1+(j2+1)*(xn1)+i+1)=1-(m*(i-x1)+y1)+j2;
%                     b2((c1-1)*(N+1)+c2)=-1;
%                     T1=[T1;b1];
%                     T2=[T2;b2];
%                     b1(1+(j1+1)*(xn1)+i+1)=0;
%                     b1(1+(j2+1)*(xn1)+i+1)=0;
%                     b2((c1-1)*(N+1)+c2)=0;
%                 end
%             end
%             if not(round(dy,5)==0)
%                 m=dx/dy;
%                 for j=round(y1):sign(y2-y1):round(y2)
%                     i1=ceil(m*(j-y1)+x1);
%                     i2=floor(m*(j-y1)+x1);
%                     b1(1+(j+1)*(xn1)+i1+1)=1-i1+(m*(j-y1)+x1);
%                     b1(1+(j+1)*(xn1)+i2+1)=1-(m*(j-y1)+x1)+i2;
%                     b2((c1-1)*(N+1)+c2)=-1;
%                     T1=[T1;b1];
%                     T2=[T2;b2];
%                     b1(1+(j+1)*(xn1)+i1+1)=0;
%                     b1(1+(j+1)*(xn1)+i2+1)=0;
%                     b2((c1-1)*(N+1)+c2)=0;
%                 end
% 
%             end
%         end
%     end
% end
% 
% 
% [s1,s3]=size(T1);
% [s4,s5]=size(T2);
% TM=sparse(s1*(tn+1),s2);
% for t=1:tn+1
%     TM((t-1)*s1+1:t*s1,(t-1)*s3+1:t*s3)=T1;
%     TM((t-1)*s1+1:t*s1,ii1+ii2+ii3+(t-1)*s5+1:ii1+ii2+ii3+t*s5)=T2;
% end
% 
end