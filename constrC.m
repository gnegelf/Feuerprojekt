function [ C ] = constrC(  inflow,capacity,N,choise)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% C=zeros(N+1,N+1,2);
% C(1,2,1)=c3;
% C(:,:,2)=ones(N+1,N+1);
% for i=1:sqrt(N)
%     for j=1:sqrt(N)
%         C((i-1)*sqrt(N)+j+1,(i-1)*sqrt(N)+j+1,1)=c1;
%         if j<sqrt(N)
%            C((i-1)*sqrt(N)+j+1,(i-1)*sqrt(N)+j+2,1)=c2;
%         end
%         if i<sqrt(N)
%            C((i-1)*sqrt(N)+j+1,(i)*sqrt(N)+j+1,1)=c2;
%         end
%     end
% end
switch choise
    case 1
        N=25;
        C=zeros(N+1,N+1,2);
        C(1,14,1)=inflow;
        c=capacity;
        %C(3,2,1)=c;
        C(8,3,1)=c;
        C(8,7,1)=c;
        C(9,8,1)=c;
        C(9,4,1)=c;
        C(14,9,1)=c;
        C(14,13,1)=c;
        C(13,12,1)=c;
        C(13,18,1)=c;
        C(18,17,1)=c;
        %C(17,22,1)=c;
        C(18,23,1)=c;
        C(14,19,1)=c;
        C(19,24,1)=c;
        C(14,15,1)=c;
        C(15,16,1)=c;
        C(15,10,1)=c;
        C(10,5,1)=c;
        C(10,11,1)=c;
        %C(11,6,1)=c;
        C(19,20,1)=c;
        C(20,21,1)=c;
        C(20,25,1)=c;
        %C(25,26,1)=c;
        C(:,:,2)=ones(N+1);
    case 2
        N=25;
        C=zeros(N+1,N+1,2);
        C(1,12,1)=inflow;
        C(1,12,2)=1;
        for i=2:sqrt(N)-2
            C(2+i,2+(i+1),1)=capacity;
            C(22+i,22+(i+1),1)=capacity;
            C(2+i,2+(i+1),2)=1;
            C(22+i,22+(i+1),2)=1;
        end

        C(12,24,1)=capacity;
        C(12,4,1)=capacity;

        C(6,11,1)=capacity;
        C(6,11,2)=1;


        C(11,16,1)=capacity;
        C(11,16,2)=1;

        C(21,16,1)=capacity;
        C(21,16,2)=1;

        C(26,21,1)=capacity;
        C(26,21,2)=1;
        C(:,:,2)=ones(N+1);
    case 3
        N=25;
        C=zeros(N+1,N+1,2);
        C(1,12,1)=inflow;
        C(1,12,2)=1;
        for i=0:sqrt(N)-2
            C(2+i,2+(i+1),1)=capacity;
            C(22+i,22+(i+1),1)=capacity;
            C(2+i,2+(i+1),2)=1;
            C(22+i,22+(i+1),2)=1;
        end
        C(7,2,1)=capacity;
        C(7,2,2)=1;
        C(6,11,1)=capacity;
        C(6,11,2)=1;

        C(12,7,1)=capacity;
        C(12,7,2)=1;
        C(11,16,1)=capacity;
        C(11,16,2)=1;
        C(12,17,1)=capacity;
        C(12,17,2)=1;
        C(21,16,1)=capacity;
        C(21,16,2)=1;
        C(17,22,1)=capacity;
        C(17,22,2)=1;
        C(26,21,1)=capacity;
        C(26,21,2)=1;
        C(:,:,2)=ones(N+1);
    case 4
        C=zeros(N+1,N+1,2);
        C(1,2,1)=inflow;
        C(:,:,2)=ones(N+1,N+1);
        for i=1:sqrt(N)
            for j=1:sqrt(N)
                %C((i-1)*sqrt(N)+j+1,(i-1)*sqrt(N)+j+1,1)=c1;
                if j<sqrt(N)
                   C((i-1)*sqrt(N)+j+1,(i-1)*sqrt(N)+j+2,1)=capacity;
                end
                if i<sqrt(N)
                   C((i-1)*sqrt(N)+j+1,(i)*sqrt(N)+j+1,1)=capacity;
                end
            end
        end
    case 5
        N=25;
        C=zeros(N+1,N+1,2);
        c=capacity;
        C(1,12,1)=inflow;
        C(12,7,1)=c;
        C(12,6,1)=c;
        C(12,13,1)=c;
        C(13,14,1)=c;
        C(13,8,1)=c;
        C(14,9,1)=c;
        C(7,8,1)=c;
        C(7,2,1)=c;
        C(2,3,1)=c;
        C(3,4,1)=c;
        C(4,5,1)=c;
        C(5,6,1)=c;
        C(8,9,1)=c;
        C(9,10,1)=c;
        C(10,11,1)=c;
        C(8,3,1)=c;
        C(9,4,1)=c;
        C(10,5,1)=c;
        C(11,6,1)=c;
        C(:,:,2)=ones(N+1,N+1);
        
end
    

end

