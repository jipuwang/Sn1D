% In this script, the approximated flux shape is tested.  It has been
% approximated with three piecewise quadratic functions. 
%% Fast flux
% Region 1:
% y = +1.5593897658139900E-03x2 -1.3533502228265900E-03x +6.3255899575709700E-02
% Region 2:
% y = -4.4443431764232400E-03x2 +2.6590845920514000E-01x -2.8376294939004300E+00
% Region 3:
% y = +1.6006807673739400E-03x2 -1.8963240555487500E-01x +5.6794357067491300E+00 

reg1 =@(x) +1.5593897658139900E-03*x.^2 ...
           -1.3533502228265900E-03*x ...
           +6.3255899575709700E-02;
reg2 =@(x) -4.4443431764232400E-03*x.^2 ...
           +2.6590845920514000E-01*x ...
           -2.8376294939004300E+00;
reg3 =@(x) +1.6006807673739400E-03*x.^2 ...
           -1.8963240555487500E-01*x ...
           +5.6794357067491300E+00;

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

%% Thermal flux
% Region 1:
% y = +3.7420664087460500E-04x2 +3.2444814149334100E-03x +1.2139249518777500E-02
% Region 2:
% y = -1.2351628643713800E-03x2 +7.3932137536755800E-02x -7.5807942857102300E-01
% Region 3:
% y = +3.6321032591887700E-04x2 -4.7314818437341200E-02x +1.5394028032516800E+00


reg1 =@(x) +3.7420664087460500E-04*x.^2 ...
           +3.2444814149334100E-03*x ...
           +1.2139249518777500E-02;
reg2 =@(x) -1.2351628643713800E-03*x.^2 ...
           +7.3932137536755800E-02*x ...
           -7.5807942857102300E-01;
reg3 =@(x) +3.6321032591887700E-04*x.^2 ...
           -4.7314818437341200E-02*x ...
           +1.5394028032516800E+00;

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