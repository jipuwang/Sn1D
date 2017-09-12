%% Instruction
  % to run different cases, change the manufacturer only!
%% Info
% Grid refiner is for grid refinement analysis. 
% Right now, it needs to be aware of the geometry and the material. So it
% can call the manufacturer to get the MMS problem and solution
% The geometry and material also need to be passed to the coupler, so the
% coupler can keep passing the info on to the modules, because it's part of
% the problem description. 
% It needs to know the geometry and is responsible for generating the grid
% and pass the grid information to the coupler. 
% function [order_phi0]=converger_SnEig(assumedSoln,assumedK)
% clear;
nGrids=4%8%4%4%6;%10;%8;
refinementRatio=2;
N=8; % angular discretization, fixed not refined. 
FDM=3;
% Geometry
Tau=10; 

% Case configure options
if ~exist('assumedSoln','var')
  assumedSoln='constant';
  assumedSoln='linear';
  assumedSoln='quadratic';
%   assumedSoln='plus1Sqrt';
%   assumedSoln='flat_expMu';
  assumedSoln='sine';
end

if ~exist('k_MMS','var')
  assumedK=1.02;
end

gridMeshSize_iGrid=zeros(nGrids,1);
error_phi0_iGrid=zeros(nGrids,1);
error_k_iGrid=zeros(nGrids,1);

for iGrid=1:nGrids
  J=5*refinementRatio^iGrid;
%   J=20;
  gridMeshSize_iGrid(iGrid)=Tau/J;
  iGrid
  % Material
  field1='Sig_t_j';          value1=ones(J,1);
  field2='Sig_ss_j';         value2=ones(J,1)*0.4;
  field3='Sig_gamma_j';      value3=ones(J,1)*0.5;
  field4='Sig_f_j';          value4=ones(J,1)*0.1;
  field5='nuSig_f_j';        value5=ones(J,1)*0.2;
  field6='thermal_cond_k_j'; value6=ones(J,1);
  field7='kappaSig_f_j';     value7=ones(J,1)*0.1; % kappa=1.0;
  mat = struct(field1,value1,field2,value2,field3,value3,... 
    field4,value4,field5,value5,field6,value6,field7,value7);

  [phi0_j_ana,psi_b1_n,psi_b2_n,Q_MMS_j_n,error_ang_j,phi0_guess_j,k_guess]=... 
        manufacturer_SnEig(J,N,Tau,mat,assumedSoln,assumedK);

  %%
%   error_ang_j=error_ang_j.*0.0;
%   Q_MMS_j_n=Q_MMS_j_n*0.0;
%   k_guess=1.0;
%   phi0_guess_j=ones(J,1);
  %%

  % Call eigen solver
  [phi0_j,k]=SnEig_module(FDM,J,N,Tau,mat,...
    psi_b1_n,psi_b2_n,Q_MMS_j_n,error_ang_j,phi0_guess_j,k_guess);

  error_phi0_iGrid(iGrid)=norm(phi0_j-phi0_j_ana-error_ang_j,2)/sqrt(J);
  error_k_iGrid(iGrid)=k-assumedK;
  
end

% Calculate the order of accuracy
order_phi0_nMinus1=ones(nGrids-1,1);
order_k_nMinus1=ones(nGrids-1,1);
for j=1:nGrids-1
  order_phi0_nMinus1(j)=log(error_phi0_iGrid(j)/error_phi0_iGrid(j+1)) / ...
    log(gridMeshSize_iGrid(j)/gridMeshSize_iGrid(j+1));
  order_k_nMinus1(j)=log(error_k_iGrid(j)/error_k_iGrid(j+1)) / ...
    log(gridMeshSize_iGrid(j)/gridMeshSize_iGrid(j+1));
end

%% Visualize the asymptotic convergence
orderPlotGrid=[gridMeshSize_iGrid(1) gridMeshSize_iGrid(end)];

scalarFluxErrorRMS_plot_handle=figure(13);
loglog(gridMeshSize_iGrid,error_phi0_iGrid,'*');
% title('scalar flux error convergence');
xlabel('mesh size [cm]');
ylabel('scalar flux error RMS');

hold on;
orderGuess=round(order_phi0_nMinus1(end));
errorStt=error_phi0_iGrid(end)*refinementRatio^(orderGuess*(nGrids-1));
firstOrder=[errorStt errorStt/refinementRatio^(nGrids-1)];
secondOrder=[errorStt errorStt/refinementRatio^(2*(nGrids-1))];
thirdOrder=[errorStt errorStt/refinementRatio^(3*(nGrids-1))];
fourthOrder=[errorStt errorStt/refinementRatio^(4*(nGrids-1))];
loglog(orderPlotGrid,firstOrder,'--');
loglog(orderPlotGrid,secondOrder,'--');
loglog(orderPlotGrid,thirdOrder,'--');
loglog(orderPlotGrid,fourthOrder,'--');
legend('scalar flux error','1st Order','2nd Order',...
  '3rd Order','4th Order','location','best');
hold off;

kErrorRMS_plot_handle=figure(14);
loglog(gridMeshSize_iGrid,error_k_iGrid,'*');
% title('scalar flux error convergence');
xlabel('mesh size [cm]');
ylabel('k error');

hold on;
orderGuess=round(order_k_nMinus1(end));
errorStt=error_phi0_iGrid(end)*refinementRatio^(orderGuess*(nGrids-1));
firstOrder=[errorStt errorStt/refinementRatio^(nGrids-1)];
secondOrder=[errorStt errorStt/refinementRatio^(2*(nGrids-1))];
thirdOrder=[errorStt errorStt/refinementRatio^(3*(nGrids-1))];
fourthOrder=[errorStt errorStt/refinementRatio^(4*(nGrids-1))];
loglog(orderPlotGrid,firstOrder,'--');
loglog(orderPlotGrid,secondOrder,'--');
loglog(orderPlotGrid,thirdOrder,'--');
loglog(orderPlotGrid,fourthOrder,'--');
legend('k error','1st Order','2nd Order',...
  '3rd Order','4th Order','location','best');
hold off;

%% Dispaly the result
% Display the problem description and results
disp '=================';
display(['assumedSoln: ' assumedSoln]);
display(['refinementRatio: ' num2str(refinementRatio)]);
display(['quad set order: ' num2str(N)]);
error_phi0_iGrid
order_phi0_nMinus1
display(num2str(order_phi0_nMinus1(end)));
error_k_iGrid
order_k_nMinus1

order_phi0=order_phi0_nMinus1(end);
order_k=order_k_nMinus1(end);

% end