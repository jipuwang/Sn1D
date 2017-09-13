% 1D MoC Module
% Input: 
%   Geometry Tau
%   Spatial discretization J (or mesh size)
%   Angular discretization N
%   Material: all cross sections
%   Boundary condition
%   Distributed source, can be MMS
% Output: 
%   Cell-averaged scalar flux

function [phi0_j,k]=SnEig_module(FDM,J,N,Tau,mat,...
           psi_b1_n,psi_b2_n,Q_MMS_j_n,...
           error_ang_j,phi0_guess_j,k_guess)
%   Input parameter
  if ~exist('FDM','var')
    FDM=3;
  end
  if ~exist('J','var')
    J=5*2;%*2%*2*2*2*2*2*2*2*2
  end
  if ~exist('N','var')
    N=16;
  end
  if ~exist('Tau','var')
    Tau=10;
  end
  if ~exist('mat','var')
    % Material
    field1='Sig_t_j';          value1=ones(J,1);
    field2='Sig_ss_j';         value2=ones(J,1)*0.5;
    field3='Sig_gamma_j';      value3=ones(J,1)*0.4;
    field4='Sig_f_j';          value4=ones(J,1)*0.1;
    field5='nuSig_f_j';        value5=ones(J,1)*0.2;
    field6='thermal_cond_k_j'; value6=ones(J,1);
    field7='kappaSig_f_j';     value7=ones(J,1)*0.1; % kappa=1.0;
    mat = struct(field1,value1,field2,value2,field3,value3,... 
      field4,value4,field5,value5,field6,value6,field7,value7);
  end
  if ~exist('psi_b1_n','var')
    psi_b1_n=ones(N,1)*0.0;
  end
  if ~exist('psi_b2_n','var')
    psi_b2_n=ones(N,1)*0.0;
  end
%   assert(norm(psi_b1_n)<1e-10,'Left incident flux is not zero.');
%   assert(norm(psi_b2_n)<1e-10,'Right incident flux is not zero.');
  if ~exist('Q_MMS_j_n','var')
    Q_MMS_j_n=ones(J,N)*0.3; % removed *2.0 (angular quantity)
  end
  if ~exist('phi0_guess_j','var')
    phi0_guess_j=ones(J,1); % For MMS problem, this cannot be arbitray.
  end
  if ~exist('k_guess','var')
    k_guess=1.0; % For MMS problem, this cannot be arbitrary.
  end
  
  % Material
  Sig_ss_j=mat.Sig_ss_j;
  nuSig_f_j=mat.nuSig_f_j;
  Sig_t_j=mat.Sig_t_j;
  Sig_t_inv_j=1./Sig_t_j;

  % Default variables, can be customized. 
  maxIterate=2000;
  outerMax=5000;
  epsilon_phi0=1e-13;
  epsilon_k=1e-10;
  delta=1E-13;
  [mu_n,weight_n]=lgwt(N,-1,1); mu_n=flipud(mu_n);
  h_j=ones(J,1)*Tau/J;
  h=Tau/J;

  % Spatial discretization Method
  switch FDM
    case 1 % Step: Step Method; Upwind; Implicit;Step; 
      alpha_n=ones(N,1);
      alpha_n(1:N/2)=-1;
    case 2 % DD: Diamond Differencing;
      alpha_n=zeros(N,1);
    case 3 % SC: Step Characteristic;
      tao_n=h./mu_n;
      alpha_n=coth(0.5*tao_n)-2./tao_n;
  end % switch FDM

  % Initial guess
  phi0_old_outer_j=phi0_guess_j;
  k_old=k_guess;
  % Initialization
  isOuterConverged=false;
  F_j=zeros(J,1);   % to store fixed fission source
  FixedSrc_j_n=zeros(J,N); % fission source plus MMS source
  
  for iOuter=1:outerMax
    % set up fission source
    F_j=1/k_old*(nuSig_f_j.*(phi0_old_outer_j));
    % calculate generalized fixed source
    for n=1:N
      FixedSrc_j_n(:,n)=F_j*0.5+Q_MMS_j_n(:,n);
    end
    % call fixed source solver
    [phi0_new_outer_j]=Sn_module(FDM,J,N,Tau,mat,...
      psi_b1_n,psi_b2_n,FixedSrc_j_n,error_ang_j,phi0_old_outer_j);
    % Be really careful here!!! The initial guess is going to be stripped
    % of the angular error in the Sn_module.  You should not do it here. 
    % update eigenvalue
    k_new=k_old*(sum(nuSig_f_j.*phi0_new_outer_j)*h)/(sum(nuSig_f_j.*phi0_old_outer_j)*h);
    %%
%     figure(111);hold on;
%     plot(phi0_new_outer_j);
%     drawnow;
    
    %%
    % check convergence
    if max(abs(phi0_new_outer_j-phi0_old_outer_j)./ ...
        (phi0_new_outer_j+delta)) < epsilon_phi0 ...
            && (abs((k_new-k_old)/k_old) <= epsilon_k)
            isOuterConverged=true;
            break;
    end
    k_old=k_new;
    phi0_old_outer_j=phi0_new_outer_j;
  end
  if ~isOuterConverged
    display (['Not converging! ' num2str(iOuter) ' outer iterations performed.']);
  else
    display (['Eigensolver converged! ' num2str(iOuter) ' outer iterations performed.']);
    phi0_j=phi0_new_outer_j;
    k=k_new;
  end
end
  
