legendList=strings(0);
timeArray=zeros(length(10:5:45)*length(30:10:60),2);
timeArray(:,1)=reshape(durations,[length(10:5:45)*length(30:10:60),1]);
durations=20000*ones(length(10:5:45),length(30:10:60));
i=1;
j=1;
for xn=10:5:10
    j=1;
    if xn==10
        for tn=30:10:60
                load(sprintf('stateNoElixn%dtn%d.mat',xn,tn));
                durations(i,j)=duration;      
            j=j+1;
        end
        i=i+1;
    else
        for tn=30:10:30
                load(sprintf('stateNoElixn%dtn%d.mat',xn,tn));
                durations(i,j)=duration;      
            j=j+1;
        end
        i=i+1;
    end
end
timeArray(:,2)=reshape(durations,[length(10:5:45)*length(30:10:60),1]);


figure
for i=1:4
    plot(10:5:45,durations(:,i),'LineWidth',3)
    set(gca,'FontSize',15)
    axis([10 45 0 20000]);
    xlabel('Spacial discritization coursety');
    ylabel('Computation time in s');
    legendList=[legendList; sprintf('Finite Differences, tn=%d',times(i))];
    legend(legendList);
    hold on;
end
savefig('noEliFigure');
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