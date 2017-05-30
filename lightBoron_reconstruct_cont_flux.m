% In this script, the approximated fast and thermal flux shape is 
% tested.  It has been approximated with three piecewise quartic
% functions. 
% y=+reg1x4*x*x*x*x +reg1x3*x*x*x +reg1x2*x*x +reg1x1*x +reg1x0;

%% For fast group
% Region 1: 
reg1x4=+1.53722295506281E-5;
reg1x3=-4.68683704257329E-4;
reg1x2=+5.77565529274875E-3;
reg1x1=-2.11263537121232E-2;
reg1x0=+2.57684318251159E-2;
% Region 2:
reg2x4=-1.10678482723517E-5;
reg2x3=+1.32004096220855E-3;
reg2x2=-6.32550721990960E-2;
reg2x1=+1.42684187372408E+0;
reg2x0=-1.14006900002113E+1;
% Region 3:
reg3x4=+3.31115176820729E-5;
reg3x3=-6.69044052900192E-3;
reg3x2=+5.06568007325930E-1;
reg3x1=-1.70499965190978E+1;
reg3x0=+2.15527557066765E+2;

reg1 =@(x) +reg1x4*x.^4 ...
           +reg1x3*x.^3 ...
           +reg1x2*x.^2 ...
           +reg1x1*x ...
           +reg1x0;
reg2 =@(x) +reg2x4*x.^4 ...
           +reg2x3*x.^3 ...
           +reg2x2*x.^2 ...
           +reg2x1*x ...
           +reg2x0;
reg3 =@(x) +reg3x4*x.^4 ...
           +reg3x3*x.^3 ...
           +reg3x2*x.^2 ...
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

% one figure to plot both
figure(101);
hold on 
fplot(reg1,[ 0 20]);
fplot(reg2,[20 40]);
fplot(reg3,[40 60]);

%% For thermal group
% Region 1: 
reg1x4=-1.75283705341027E-5;
reg1x3=+5.77040031150232E-4;
reg1x2=-5.12351183707553E-3;
reg1x1=+3.64338769502884E-2;
reg1x0=-5.53011823387657E-3;
% Region 2:
reg2x4=+1.98913688295027E-5;
reg2x3=-2.38197650467746E-3;
reg2x2=+1.05660411771455E-1;
reg2x1=-2.05691181513866E+0;
reg2x0=+1.52327882685776E+1;
% Region 3:
reg3x4=-1.86512366899547E-5;
reg3x3=+3.86918856164260E-3;
reg3x2=-2.98883311241798E-1;
reg3x1=+1.01579055107311E+1;
reg3x0=-1.27495442059681E+2;

reg1 =@(x) +reg1x4*x.^4 ...
           +reg1x3*x.^3 ...
           +reg1x2*x.^2 ...
           +reg1x1*x ...
           +reg1x0;
reg2 =@(x) +reg2x4*x.^4 ...
           +reg2x3*x.^3 ...
           +reg2x2*x.^2 ...
           +reg2x1*x ...
           +reg2x0;
reg3 =@(x) +reg3x4*x.^4 ...
           +reg3x3*x.^3 ...
           +reg3x2*x.^2 ...
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

% one figure to plot both
figure(101);
fplot(reg1,[ 0 20]);
fplot(reg2,[20 40]);
fplot(reg3,[40 60]);
title('reconstructed fast and thermal flux');
legend('fast region 1','fast region 2','fast region 3',...
  'thermal region 1','thermal region 2','thermal region 3');
hold off;
grid on;
