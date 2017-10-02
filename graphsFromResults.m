legendList=strings(0);
timelimit=20000;
scenario=2;


cbStruct.namecode='Results/state';
cbStruct.xn=10:5:15;
cbStruct.tn=30:10:40;
cbStruct.durations=timelimit*ones(length(cbStruct.xn),length(cbStruct.tn));
cbStruct.cbnums=zeros(length(cbStruct.xn),length(cbStruct.tn));
cbStruct.constrnums=zeros(length(cbStruct.xn),length(cbStruct.tn));
cbStruct.legendLabel='Finite Elements w Callback';
cbStruct.figName='FEwCB';
cbStruct.axis=[min(cbStruct.xn),max(cbStruct.xn),1,20000];


noCbStruct.namecode='Results/stateFull';
noCbStruct.xn=10:5:15;
noCbStruct.tn=30:10:40;
noCbStruct.durations=timelimit*ones(length(noCbStruct.xn),length(noCbStruct.tn));
noCbStruct.cbnums=zeros(length(noCbStruct.xn),length(noCbStruct.tn));
noCbStruct.constrnums=zeros(length(noCbStruct.xn),length(noCbStruct.tn));
noCbStruct.legendLabel='Finite Elements w/o Callback';
noCbStruct.figName='FEwoCB';
noCbStruct.axis=[min(noCbStruct.xn),max(noCbStruct.xn),1,20000];


% noEliStruct.namecode='Results/stateNoEli';
% noEliStruct.xn=10:5:10;
% noEliStruct.tn=30:10:40;
% noEliStruct.durations=timelimit*ones(length(noEliStruct.xn),length(noEliStruct.tn));
% noEliStruct.cbnums=zeros(length(noEliStruct.xn),length(noEliStruct.tn));
% noEliStruct.constrnums=zeros(length(noEliStruct.xn),length(noEliStruct.tn));
% noEliStruct.legendLabel='Finite Differences';
% noEliStruct.figName='FDnoEli';
% noEliStruct.axis=[min(noEliStruct.xn),max(noEliStruct.xn),1,20000];

%cbStruct.xlabel='Spacial discretization resolution';
%cbStruct.ylabel='Computation time in s';
models={noCbStruct,cbStruct};

for k=1:length(models)
    legendList=strings(0);
    models{k}.xlabel='Spacial discretization resolution';
    models{k}.ylabel='Computation time in s';
    for i=1:length(models{k}.xn)
        for j=1:length(models{k}.tn)
            data=load(sprintf([models{k}.namecode 'xn%dtn%ds%d.mat'],models{k}.xn(i),models{k}.tn(j),scenario));
            if isfield(data,'duration')
                models{k}.durations(i,j)=data.duration;
            end
            if isfield(data,'cbnum')
                models{k}.cbnums(i,j)=data.cbnum;
            end
            if isfield(data,'number_of_constraints')
               models{k}.constrnums(i,j)=data.number_of_constraints; 
            end
        end

    end

    fig=figure
    for t=1:length(models{k}.tn)
        semilogy(models{k}.xn,models{k}.durations(:,t),'LineWidth',3)
        set(gca,'FontSize',18)
        axis(models{k}.axis);
        xlabel(models{k}.xlabel,'FontSize',22);
        ylabel(models{k}.ylabel,'FontSize',22);
        legendList=[legendList; sprintf([models{k}.legendLabel,' tn=%d'],models{k}.tn(t))];
        hold on;
    end
    legend(legendList,'FontSize',17,'Location','southeast');
    %saveas(fig,'figureFiniteDifferences.eps');
    print(models{k}.figName,'-depsc ');
    
end
save('plotData','models');
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