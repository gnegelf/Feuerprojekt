function [ C,contPos,N ] = constrGraph(  inflow,capacity,choice)
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
switch choice
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
        contPos=zeros(N,2);
        for i=1:sqrt(N)
            for j=1:sqrt(N)
                contPos((i-1)*sqrt(N)+j,1)=1/(sqrt(N)+1)*j;
                contPos((i-1)*sqrt(N)+j,2)=1/(sqrt(N)+1)*i;
            end
        end
    case 6
        N=14;
        C=zeros(N+1,N+1,2);
        c=capacity;
        C(1,2,1)=inflow;
        C(2,[3,4,5,6],1)=c;
        C([3,9],8,1)=c;
        C([4,8,10],9,1)=c;
        C([5,9,11],10,1)=c;
        C(6,7,1)=c;
        C([7,10],11,1)=c;
        C(8,12,1)=c;
        C(9,13,1)=c;
        C(10,14,1)=c;
        C(11,15,1)=c;
        C(:,:,2)=ones(N+1,N+1);
        contPos=[0.1,0.5;0.2,0.9;0.2,0.6;0.2,0.4;0.2,0.1;0.4,0.1;0.75,0.9;0.75,0.6;0.75,0.4;0.75,0.1;0.9,0.9;0.9,0.6;0.9,0.4;0.9,0.1];
    case 7
        N=7;
        C=zeros(N+1,N+1,2);
        c=capacity;
        C(1,2,1)=inflow;
        C(2,[3,8],1)=c;
        C(3,4,1)=c;
        C(4,5,1)=c;
        C(5,6,1)=c;
        C([6,8],7,1)=c;
        C(2,8,1)=c;
        C(:,:,2)=ones(N+1,N+1);
        contPos=[0.15,0.15;0.15,0.5;0.15,0.85;0.4,0.85;0.6,0.85;0.85,0.85;0.85,0.15];
    case 8
        N=18;
        C=zeros(N+1,N+1,2);
        c=capacity;
        C(1,2,1)=inflow;
        C(2,[3,15],1)=c;
        C(3,4,1)=c;
        C(4,5,1)=c;
        C(5,6,1)=c;
        C(6,7,1)=c;
        C(7,8,1)=c;
        C(8,9,1)=c;
        C(9,10,1)=c;
        C(10,11,1)=c;
        C(11,12,1)=c;
        C(12,13,1)=c;
        C(13,14,1)=c;
        C(15,[16,18],1)=c;
        C(16,17,1)=c;
        C(18,19,1)=c;
        C(:,:,2)=ones(N+1,N+1);
        contPos=[0.05,0.05; 0.05,0.9; 0.3,0.9; 0.4,0.9; 0.5,0.9; 0.55,0.9; 0.6,0.9; 0.65,0.9; 0.7,0.9; 0.8,0.9; 0.9,0.9; 0.9,0.75; 0.9,0.6; 0.4,0.4; 0.45,0.6; 0.45,0.8; 0.9,0.4; 0.95 ,0.7];
    case 9
        N=17;
        C=zeros(N+1,N+1,2);
        c=capacity;
        C(1,2,1)=inflow;
        C(2,[3,7,11,15],1)=c;
        for i=[3,4,5,7,8,9,11,12,13,15,16,17]
           C(i,i+1,1)=c; 
        end
        C(:,:,2)=ones(N+1,N+1);
        contPos=[0.2,0.5; 0.6,0.2; 0.7,0.2; 0.8,0.2; 0.9,0.2; 0.6,0.4; 0.7,0.4; 0.8,0.4; 0.9,0.4; 0.6,0.6; 0.7,0.6; 0.8,0.6; 0.9,0.6; 0.6,0.8; 0.7,0.8; 0.8,0.8; 0.9,0.8];
    case 10
        N=18;
        C=zeros(N+1,N+1,2);
        c=capacity;
        C(1,2,1)=inflow;
        C(1,19,1)=inflow;
        C(2,[3,7,11,15],1)=c;
        C(19,[6,10,14,18],1)=c;
        for i=[3,4,5,7,8,9,11,12,13,15,16,17]
           C(i,i+1,1)=c; 
        end
        for i=[4,5,8,9,12,13,16,17,6,10,14,18]
           C(i,i-1,1)=c; 
        end
        C(:,:,2)=ones(N+1,N+1);
        contPos=[0.1,0.5; 0.25,0.2; 0.4,0.2; 0.6,0.2; 0.75,0.2; 0.25,0.4; 0.4,0.4; 0.6,0.4; 0.75,0.4; 0.25,0.6; 0.4,0.6; 0.6,0.6; 0.75,0.6; 0.25,0.8; 0.4,0.8; 0.6,0.8; 0.75,0.8; 0.9,0.5];
    case 11
        N=36;
        Ns=sqrt(N);
        C=zeros(N+1,N+1,2);
        contPos=zeros(N,2);
        for i=1:Ns
            for j=1:Ns
                contPos((i-1)*Ns+j,1)=i/(Ns+1);
                contPos((i-1)*Ns+j,2)=j/(Ns+1);
            end
        end
    case 12
        N=36;
        C=zeros(N+1,N+1,2);
        contPos=zeros(N,2);
        for i=1:N
                contPos(i,1)=0.1;
                contPos(i,2)=i/(N+1);
        end
end
    
        
    

end

