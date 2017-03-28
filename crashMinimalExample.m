figure
plot([0,1],[0,1]);
hold on;
plot([0,1],[0.5,0.5]);
hold on;
axis([0,1,0,1]);
xlabel('x');
ylabel('y');
legendList=strings(2);
legendList(1)='1';
legendList(2)='2';
legend(legendList);