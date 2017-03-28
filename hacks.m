times=30:10:60;
spaces=10:5:45;
duration=zeros(length(10:5:45),length(30:10:60));
i=1;
j=1;
legendList=strings(0);
for xn=10:5:45
    j=1;
    for tn=30:10:60
        load(sprintf('statexn%dtn%d.mat',xn,tn));
        duration(i,j)=eval(sprintf('duration%d',xn));
        j=j+1;
    end
    i=i+1;
end

figure
for i=1:4
    pp=plot(10:5:45,duration(:,i));
    axis([10 45 0 20000]);
    xlabel('Spacial discritization coursety');
    ylabel('Computation time in s');
    legendList=[legendList; sprintf('elimination, tn=%d',times(i))];
    %clabel('Hallo');
    hold on;
end
legend(legendList);
legendList=strings(0);

durations=20000*ones(length(10:5:45),length(30:10:60));
i=1;
j=1;
for xn=10:5:45
    j=1;
    for tn=30:10:60
        if (tn<45 || (tn==50 && xn<42) || (tn==60 && xn < 31))
            load(sprintf('stateFullxn%dtn%d.mat',xn,tn));
            durations(i,j)=duration;
        end
        
        j=j+1;
    end
    i=i+1;
end



figure
for i=1:4
    plot(10:5:45,durations(:,i))
    axis([10 45 0 20000]);
    xlabel('Spacial discritization coursety');
    ylabel('Computation time in s');
    legendList=[legendList; sprintf('full, tn=%d',times(i))];
    legend(legendList);
    hold on;
end
clear variables;
% 
% figure
% plot([0,1],[0,1]);
% hold on;
% plot([0,1],[0.5,0.5]);
% hold on;
% axis([0,1,0,1]);
% xlabel('x');
% ylabel('y');
% legendList=strings(2);
% legendList(1)='1';
% legendList(2)='2';
% legend(legendList);