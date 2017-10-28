slow=2;
if eliminate
    z=solBasis*x_k(1:i2)+solInhom;
    plotM=reshape(z,[xn1,xn1,tn1]);
    plotMEli=plotM;
else
    plotM=reshape(x_k(1:(xn+3)^2*tn1),[xn+3,xn+3,tn1]);
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
p=x_k(1:1+tn*N);%p enth�lt Eintr�ge der Wasserentnahme
extractMatrix=zeros(N,tn);
for t=1:tn
    for i=1:N
        extractMatrix(i,t)=p((t-1)*N+i); 
    end
end

fig=figure;
plot([0 1], [0 0],'b',[0 1], [1 1],'b')
axis tight manual
ax = gca;
ax.NextPlot = 'replaceChildren';

%rescale temperature of the tiles in plot to the mean of the four cornors

loops = tn+1;
clear F;
F(loops*slow) = struct('cdata',[],'colormap',[]);
clear vid;
%q=x_k(i2+2*i3+1:length(x_k));%q enth�lt Ein
vid((loops))=struct('cdata',[],'colormap',[]);
sharp=0;
for j = 1:loops
    for jj=1:slow
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
        plot(Graph(:,1),Graph(:,2),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'MarkerSize',4);
        for i=1:N
            if j > 1
                val=0;
                val=val+extractMatrix(i,j-1);
                if val>0
                    plot(Graph(i,1),Graph(i,2),'d','MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'MarkerSize',val*10);%wasserausguss
                end
            end
        end
        drawnow
        vid((j-1)*slow+jj) = getframe(fig);
       F((j-1)*slow+jj)=getframe(fig);   
    end
end
%movie(G,1,40);
str1=sprintf('testEli:%d-x:%d-t:%d',eliminate,xn,tn);
savevar=0;
if savevar
v=VideoWriter(str1);
open(v);
writeVideo(v,F);
close(v);
end
end