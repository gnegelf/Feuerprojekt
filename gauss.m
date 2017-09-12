gausss=@(x,y) exp(-50*((x-0.5)^2+(y-0.5)^2));
x=0:0.025:1;
plotM=zeros(length(x),length(x));
for i=1:length(x)
    for j=1:length(x)
        plotM(i,j)=gausss((i-1)/length(x),(j-1)/length(x));
    end
end
fig=figure;
surf(x,x,0.1*plotM);
axis([0 1 0 1 0 1]);
print('gauss01','-depsc');
fig1 = figure;

surf(x,x,0.5*plotM);
axis([0 1 0 1 0 1]);
print('gauss05','-depsc');
fig2=figure;
surf(x,x,1*plotM);
axis([0 1 0 1 0 1]);
print('gauss1','-depsc');
