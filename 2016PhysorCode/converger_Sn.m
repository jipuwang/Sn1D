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
function [order_phi]=converger_Sn(assumedSoln,FDM)
% clear;
nGrids=4%8%4%4%6;%10;%8;
refinementRatio=2;
N=2; % angular discretization, fixed not refined. 

% Geometry
Tau=10; 

% Case configure options
if ~exist('FDM','var')
  FDM=3;
end
if ~exist('assumedSoln','var')
  assumedSoln='constant';
  assumedSoln='linear';
  assumedSoln='quadratic';
%   assumedSoln='plus1Sqrt';
%   assumedSoln='flat_expMu';
end

error_phi0_n=zeros(nGrids,1);
gridMeshSize_n=zeros(nGrids,1);

for iGrid=1:nGrids
  J=5*refinementRatio^iGrid;
  gridMeshSize_n(iGrid)=Tau/J;
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

  [phi0_j_ana,psi_b1_n,psi_b2_n,Q_MMS_j_n,error_ang_j]=... 
        manufacturer_Sn(J,N,Tau,mat,assumedSoln);
%   error_ang_j=error_ang_j*0.0;
  [phi0_j]=Sn_module(FDM,J,N,Tau,mat,...
    psi_b1_n,psi_b2_n,Q_MMS_j_n,error_ang_j);

  % Calculate the error compared to manufactured solution
%   error_ang_j=zeros(J,1);
% error_ang_j=error_ang_j.*0.0;
  error_phi0_n(iGrid)=norm(phi0_j-phi0_j_ana-error_ang_j,2)/sqrt(J) 
  
end

% Calculate the order of accuracy
order_phi_nMinus1=ones(nGrids-1,1);
for j=1:nGrids-1
  order_phi_nMinus1(j)=log(error_phi0_n(j)/error_phi0_n(j+1)) / ...
    log(gridMeshSize_n(j)/gridMeshSize_n(j+1));
end

%% Visualize the asymptotic convergence
orderPlotGrid=[gridMeshSize_n(1) gridMeshSize_n(end)];

scalarFluxErrorRMS_plot_handle=figure(13);
loglog(gridMeshSize_n,error_phi0_n,'*');
title({'cell-averaged scalar flux error convergence',...
  ['\phi_{MMS}: ' assumedSoln]});
xlabel('mesh size [cm]');
ylabel('scalar flux error RMS');

hold on;
orderGuess=round(order_phi_nMinus1(end));
errorStt=error_phi0_n(end)*refinementRatio^(orderGuess*(nGrids-1));
firstOrder=[errorStt errorStt/refinementRatio^(nGrids-1)];
secondOrder=[errorStt errorStt/refinementRatio^(2*(nGrids-1))];
thirdOrder=[errorStt errorStt/refinementRatio^(3*(nGrids-1))];
fourthOrder=[errorStt errorStt/refinementRatio^(4*(nGrids-1))];
loglog(orderPlotGrid,firstOrder,'--');
loglog(orderPlotGrid,secondOrder,'--');
loglog(orderPlotGrid,thirdOrder,'--');
loglog(orderPlotGrid,fourthOrder,'--');
legend('scalar flux error','1st Order','2nd Order',...
  '3rd Order','4th Order','location','northwest');
hold off;
% savefig(scalarFluxErrorRMS_plot_handle,['temp_Sn_cellAveraged_phi_convergence_' assumedSoln]);
savefig(scalarFluxErrorRMS_plot_handle,['Sn_Physor2016_cellAveraged_phi_convergence_' assumedSoln]);
%% Dispaly the result
% Display the problem description and results
disp '=================';
display(['assumedSoln: ' assumedSoln]);
display(['Number of grids: ' num2str(nGrids)]);
display(['refinementRatio: ' num2str(refinementRatio)]);
display(['FDM method: ' num2str(FDM)]);
display(['quad set order: ' num2str(N)]);
error_phi0_n
order_phi_nMinus1
display(num2str(order_phi_nMinus1(end)));
order_phi=order_phi_nMinus1(end);

end