slow=5;


%z=solInhom;
if ~finiteDifferences
    z=solBasis*x_k(1:i2)+solInhom;
    plotM=reshape(z,[xn1,xn1,tn1]);
    plotMEli=plotM;
else
    plotM=reshape(Result.x_k(1:(xn+3)^2*tn1),[xn+3,xn+3,tn1]);
    plotM=plotM(2:xn+2,2:xn+2,:);
    plotMNoEli=plotM;
end

fplotM=zeros(size(plotM));
for t=1:tn+1
    for i=1:xn+1
        for j=1:xn+1
            if i<xn+1 && j<xn+1
                fplotM(i,j,t)=(plotM(i,j,t)+plotM(i+1,j,t)+plotM(i,j+1,t)+plotM(i+1,j+1,t))/4;
            else
                fplotM(i,j,t)=plotM(i,j,t);
            end
        end
    end
end
if video
x=0:dx:1;
y=x;
Graph=G;
fine=8;
flowMatrix=zeros(N,N,tn+1);
p=x_k(tn*N+1:length(x_k));% p enth�lt die Eintr�ge der Flussvariablen
for t=1:tn+1
    for i=1:N
        for j=1:N
            flowMatrix(i,j,t)=p((t-1)*(N+1)^2+i*(N+1)+j+1);
        end
    end
end
%flowMatrix enth�lt die Flussvariablen in Matrixform
p=x_k(1:1+tn*N);%p enth�lt Eintr�ge der Wasserentnahme
extractMatrix=zeros(N,tn);
for t=1:tn
    for i=1:N
        extractMatrix(i,t)=p((t-1)*N+i); 
    end
end
%extractMatrix enth�lt die Wasserentnahme an den Knoten am jeweiligen
%Zeitschritt

fig=figure;
plot([0 1], [0 0],'b',[0 1], [1 1],'b')
axis tight manual
ax = gca;
ax.NextPlot = 'replaceChildren';
useStepAmount=ceil(usetime/dt);

%rescale temperature of the tiles in plot to the mean of the four cornors

loops = tn+1;
clear F;
F(3*slow*(loops)) = struct('cdata',[],'colormap',[]);
clear vid;
p=x_k(tn*N+1:length(x_k));%p enth�lt Eintr�ge ohne Fluss und Entnahmevariablen
%q=x_k(i2+2*i3+1:length(x_k));%q enth�lt Ein
vid(slow*(loops))=struct('cdata',[],'colormap',[]);
stredge=15/fac;
stredge2=30/fac;
stredge3=40/fac;
sharp=1;

for j = 1:loops
    for ii=1:slow
        clf
        str1=sprintf('Time: %0.2f',(j-1)*dt);
        plot([0 1], [0 0],'k',[0 1], [1 1],'k')
        title(str1);
        hold on;
        plot([0 0], [0 1],'k',[1 1], [0 1],'k')
        if sharp
            pcolor(x,y,(plotM(:,:,j)>1.0)*1.0)
            colorbar
            caxis([0,max(max(max(plotM>1.0)))])
        else
            pcolor(x,y,(plotM(:,:,j)))
            colorbar
            caxis([0,max(max(max(plotM)))])
        end
        
        for i=1:N
           for k=1:N
               if C(i+1,k+1,1)>0
                    line([Graph(i,1),Graph(k,1)],[Graph(i,2),Graph(k,2)],'Color',[1 1 1],'LineWidth',stredge*C(i+1,k+1,1));
               end
           end
        end
        for i=1:N
           for k=1:N
               for cc=0:C(i+1,k+1,2)-1               
                    for co =0:step-1
                        if j>co+cc*step
                            if flowMatrix(i,k,j-cc*step-co)>10^(-5)
                                x1=Graph(i,1)+(Graph(k,1)-Graph(i,1))*((ii+cc*slow*step+co*slow-1)/(slow*step*C(i+1,k+1,2)));
                                y1=Graph(i,2)+(Graph(k,2)-Graph(i,2))*((ii+cc*slow*step+co*slow-1)/(slow*step*C(i+1,k+1,2)));
                                x2=Graph(i,1)+(Graph(k,1)-Graph(i,1))*((ii+cc*slow*step+co*slow)/(slow*step*C(i+1,k+1,2)));
                                y2=Graph(i,2)+(Graph(k,2)-Graph(i,2))*((ii+cc*slow*step+co*slow)/(slow*step*C(i+1,k+1,2)));
                                line([x1,x2],[y1,y2],'Color',[0 0 0],'LineWidth',stredge*flowMatrix(i,k,j-cc*step-co));
                            end
                        end
                    end
               end
           end
           if p((j-1)*(N+1)^2+i+1)>0

               %plot(Graph(i,1),Graph(i,2),'s','MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 1 1],'MarkerSize',stredge2*p((j-1)*(N+1)^2+i+1));%flow from source
           end
        end
        plot(Graph(:,1),Graph(:,2),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'MarkerSize',4);
        for i=1:N
            if j > 1
                val=extractMatrix(i,j-1);
                %for k=0:useStepAmount-1
                %    if(j-1-k)>0
                %        val=val+extractMatrix(i,j-1-k);
                %    end
                %end
                if val>0
                    plot(Graph(i,1),Graph(i,2),'d','MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'MarkerSize',stredge3*val);%wasserausguss
                end
                
            end
        end
        %plot(Graph(:,1),Graph(:,2),'o','MarkerEdgeColor',[1 0 0],'MarkerSize',4);
        drawnow
        vid((j-1)*slow+ii) = getframe(fig);
    for kk=1:3
       F((j-1)*slow*3+(ii-1)*3+kk)=getframe(fig);
    end
    end
    
end
%movie(G,1,40);
str1=sprintf('testEli%d-x%d-t%d.mp4',eliminate,xn,tn);
v=VideoWriter(str1);
open(v);
writeVideo(v,F);
close(v);

totalFlow=zeros(tn,1);
for i=2:tn+1
    totalFlow(i-1)=sum(sum(flowMatrix(:,:,i)));
end

end