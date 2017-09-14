legendList=strings(0);
%timeArray=zeros(length(10:5:45)*length(30:10:60),2);
space=[10:1:15,20:5:45];
durations=20000*ones(length(space),length(30:10:60));
%timeArray(:,1)=reshape(durations,[length(10:5:45)*length(30:10:60),1]);
times=[30,40,50,60];
i=1;
j=1;
for xn=10:1:12
    j=1;
    for tn=30:10:60
            load(sprintf('stateNoElixn%dtn%d.mat',xn,tn));
            durations(i,j)=duration;      
        j=j+1;
    end
    i=i+1;
end
j=1;
for tn=[30,40,50,60]
    if tn ~= 50
            load(sprintf('stateNoElixn%dtn%d.mat',13,tn));
            durations(i,j)=duration; 
    end
        j=j+1;
end
i=i+1;
j=4;
load(sprintf('stateNoElixn%dtn%d.mat',14,60));
durations(i,j)=duration;      
%timeArray(:,2)=reshape(durations,[length(10:5:45)*length(30:10:60),1]);


fig=figure
for i=1:4
    semilogy(space,durations(:,i),'LineWidth',3)
    set(gca,'FontSize',18)
    axis([10 45 600 20000]);
    xlabel('Spacial discretization resolution','FontSize',22);
    ylabel('Computation time in s','FontSize',22);
    legendList=[legendList; sprintf('Finite Differences, tn=%d',times(i))];
    hold on;
end
axis([10 45 600 20000]);
legend(legendList,'FontSize',17,'Location','southeast');
%saveas(fig,'figureFiniteDifferences.eps');
print -depsc figureFiniteDifferences;
save('durationsNoEli','durations');
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