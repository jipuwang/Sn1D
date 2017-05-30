% In this script, the approximated both fast and thermal flux shape is 
% tested.  It has been approximated with three piecewise quadratic
% functions. 
% y=reg1x2*x*x +reg1x1*x +reg1x0;

%% For fast group
% Region 1: 
reg1x2=+0.0020912236785901;
reg1x1=-0.0088038653251593;
reg1x0=+0.0345867334422320;
% Region 2:
reg2x2=-0.0044008146496444;
reg2x1=+0.2639624489114600;
reg2x0=-2.8239242199963540;
% Region 3:
reg3x2=+0.0019059413144357;
reg3x1=-0.2235676470181600;
reg3x0=+6.5864700746602860;

reg1 =@(x) +reg1x2*x.^2 ...
           +reg1x1*x ...
           +reg1x0;
reg2 =@(x) +reg2x2*x.^2 ...
           +reg2x1*x ...
           +reg2x0;
reg3 =@(x) +reg3x2*x.^2 ...
           +reg3x1*x ...
           +reg3x0;

figure(99); 
hold on;
fplot(reg1,[ 0 20]);
fplot(reg2,[20 40]);
fplot(reg3,[40 60]);
title('reconstructed fast flux');
legend('region 1','region 2','region 3');
hold off;
grid on;

% to fit everything in the same plot
figure(101); 
hold on;
fplot(reg1,[ 0 20]);
fplot(reg2,[20 40]);
fplot(reg3,[40 60]);

%% For thermal group
% Region 1: 
reg1x2=+0.00038256921789935;
reg1x1=+0.00372575820304040;
reg1x0=+0.00245833122976550;
% Region 2:
reg2x2=-0.00116713935976200;
reg2x1=+0.07017666404691600;
reg2x0=-0.70667635458320700;
% Region 3:
reg3x2=+0.00038371079460677;
reg3x1=-0.04964502056042800;
reg3x0=+1.60483078272052200;

reg1 =@(x) +reg1x2*x.^2 ...
           +reg1x1*x ...
           +reg1x0;
reg2 =@(x) +reg2x2*x.^2 ...
           +reg2x1*x ...
           +reg2x0;
reg3 =@(x) +reg3x2*x.^2 ...
           +reg3x1*x ...
           +reg3x0;

figure(100); 
hold on;
fplot(reg1,[ 0 20]);
fplot(reg2,[20 40]);
fplot(reg3,[40 60]);
title('reconstructed thermal flux');
legend('region 1','region 2','region 3');
hold off;
grid on;

% to fit everything in the same plot
figure(101); 
fplot(reg1,[ 0 20]);
fplot(reg2,[20 40]);
fplot(reg3,[40 60]);
title('reconstructed fast and thermal flux');
legend('fast region 1','fast region 2','fast region 3',...
  'thermal region 1','thermal region 2','thermal region 3');
hold off;
grid on;