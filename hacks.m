times=30:10:60;
spaces=10:5:45;
durations=zeros(length(10:5:45),length(30:10:60));
i=1;
legendList=strings(0);
for xn=10:5:45
    j=1;
    for tn=30:10:60
        load(sprintf('statexn%dtn%d.mat',xn,tn));
        durations(i,j)=eval(sprintf('duration',xn));
        j=j+1;
    end
    i=i+1;
end

figure 
for i=1:4
    pp=semilogy(10:5:45,durations(:,i),'LineWidth',3);
    set(gca,'FontSize',18)
    axis([10 45 1 20000]);
    xlabel('Spacial discritization resolution','FontSize',22);
    ylabel('Computation time in s','FontSize',22);
    legendList=[legendList; sprintf('With Callback, tn=%d',times(i))];
    %clabel('Hallo');
    hold on;
end
save('durationsLazy','durations');
legend(legendList,'FontSize',17,'Location','southeast');
legendList=strings(0);

savefig('eliminationFigure');
print -depsc figureReducedModelCB;
%print -png figureReducedModelCB;
timeArray=zeros(length(10:5:45)*length(30:10:60),2);
timeArray(:,1)=reshape(durations,[length(10:5:45)*length(30:10:60),1]);
durations=20000*ones(length(10:5:45),length(30:10:60));
i=1;
j=1;
for xn=10:5:45
    j=1;
    for tn=30:10:60
            load(sprintf('stateFullxn%dtn%d.mat',xn,tn));
            durations(i,j)=duration;      
        j=j+1;
    end
    i=i+1;
end
timeArray(:,2)=reshape(durations,[length(10:5:45)*length(30:10:60),1]);


fig = figure
for i=1:4
    semilogy(10:5:45,durations(:,i),'LineWidth',3)
    set(gca,'FontSize',18)
    axis([10 45 1 20000]);
    xlabel('Spacial discretization resolution','FontSize',22);
    ylabel('Computation time in s','FontSize',22);
    legendList=[legendList; sprintf('Finite Elements, tn=%d',times(i))];
    
    hold on;
end
legend(legendList,'FontSize',17,'Location','southeast');

%saveas(fig,'figureReducedModel.png');
print -depsc figureReducedModel;
%print -png figureReducedModel;
save('zeiten','timeArray');
save('durationsFull','durations');
% 
% figure
% plot([0,1],[0,1]);
% hold on;
% plot([0,1],[0.5,0.5]);length(10:5:45),length(30:10:60)
% hold on;
% axis([0,1,0,1]);
% xlabel('x');
% ylabel('y');
% legendList=strings(2);
% legendList(1)='1';
% legendList(2)='2';
% legend(legendList);