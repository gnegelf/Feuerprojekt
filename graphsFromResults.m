legendList=strings(0);
timelimit=20000;
scenario=2;


cbStruct.namecode='Results/state';
cbStruct.xn=[10:5:45 50:10:70 80:20:160];
cbStruct.tn=30:10:60;
cbStruct.durations=timelimit*ones(length(cbStruct.xn),length(cbStruct.tn));
cbStruct.cbnums=zeros(length(cbStruct.xn),length(cbStruct.tn));
cbStruct.constrnums=zeros(length(cbStruct.xn),length(cbStruct.tn));
cbStruct.legendLabel='Finite Elements w Callback';
cbStruct.figName='FEwCB';
cbStruct.axis=[min(cbStruct.xn),max(cbStruct.xn),1,20000];


noCbStruct.namecode='Results/stateFull';
noCbStruct.xn=[10:5:45 50:10:80 100:20:160];
noCbStruct.tn=30:10:60;
noCbStruct.durations=timelimit*ones(length(noCbStruct.xn),length(noCbStruct.tn));
noCbStruct.cbnums=zeros(length(noCbStruct.xn),length(noCbStruct.tn));
noCbStruct.constrnums=zeros(length(noCbStruct.xn),length(noCbStruct.tn));
noCbStruct.legendLabel='Finite Elements w/o Callback';
noCbStruct.figName='FEwoCB';
noCbStruct.axis=[min(noCbStruct.xn),max(noCbStruct.xn),1,20000];


noEliStruct.namecode='Results/stateNoEli';
noEliStruct.xn=10:1:13;
noEliStruct.tn=30:10:60;
noEliStruct.durations=timelimit*ones(length(noEliStruct.xn),length(noEliStruct.tn));
noEliStruct.cbnums=zeros(length(noEliStruct.xn),length(noEliStruct.tn));
noEliStruct.constrnums=zeros(length(noEliStruct.xn),length(noEliStruct.tn));
noEliStruct.legendLabel='Finite Differences';
noEliStruct.figName='FDnoEli';
noEliStruct.axis=[min(noEliStruct.xn),max(noEliStruct.xn),1,20000];

%cbStruct.xlabel='Spacial discretization resolution';
%cbStruct.ylabel='Computation time in s';
models={noCbStruct,cbStruct,noEliStruct};

for k=1:length(models)
    legendList=strings(0);
    models{k}.xlabel='Spacial discretization resolution';
    models{k}.ylabel='Computation time in s';
    for i=1:length(models{k}.xn)
        for j=1:length(models{k}.tn)
            try
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
            catch
                
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
legendList2=strings(0);
for tn=30:10:60
   if tn==30
       t=1
   end
   if tn==40
       t=2
   end
   if tn==50
       t=3
   end
   if tn==60
       t=4
   end
           
           
   legendList2=strings(0);
   fig2=figure
   semilogy(models{1}.xn,models{1}.durations(:,t),'LineWidth',3)
   hold on
   semilogy(models{2}.xn,models{2}.durations(:,t),'LineWidth',3)
   hold on
   semilogy(models{3}.xn,models{3}.durations(:,t),'LineWidth',3)
   title(sprintf('%d timesteps',tn))
   set(gca,'FontSize',18)
   axis([10,45,1,20000]);
   xlabel(models{1}.xlabel,'FontSize',22);
   ylabel(models{1}.ylabel,'FontSize',22);
   legendList2=[legendList2; models{1}.legendLabel;models{2}.legendLabel;models{3}.legendLabel];
   
   hold off;
   legend(legendList2,'FontSize',17,'Location','southeast');
end
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