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

function [phi0_j]=Sn_module(J,N,Tau,mat,...
           psi_b1_n,psi_b2_n,Q_MMS_j_n,error_ang_j)

%   Input parameter
  if ~exist('Tau','var')
    Tau=10;
  end
  if ~exist('J','var')
    J=5*2;%*2%*2*2*2*2*2*2*2*2
  end
  if ~exist('N','var')
    N=16;
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
    psi_b1_n=ones(N,1)*1.0;
  end
  if ~exist('psi_b2_n','var')
    psi_b2_n=ones(N,1)*1.0;
  end
  if ~exist('Q_MMS_j_n','var')
    Q_MMS_j_n=ones(J,N)*0.3; % removed *2.0 (angular quantity)
  end
  
  % Material
  Sig_ss_j=mat.Sig_ss_j;
  nuSig_f_j=mat.nuSig_f_j;
  Sig_t_j=mat.Sig_t_j;

  Sig_t_inv_j=1./Sig_t_j;
  
  % Default variables, can be customized. 
  maxIterate=2000;
  FDM=3;
  epsilon_phi0=1e-13;
  delta=1E-13;
  [mu_n,weight_n]=lgwt(N,-1,1); mu_n=flipud(mu_n);
  h_j=ones(J,1)*Tau/J;
  h=Tau/J;
  
  % Spatial discretization Method
  % =1 -> Step Method; Upwind; Implicit;Step;
  % =2 -> Diamond Differencing;DD;
  % =3 -> Step Characteristic;SC;
  % would be alpha_nj(N,J) with nonuniform spatial grid
  switch FDM
    case 1 % Step
      alpha_n=ones(N,1);
      alpha_n(1:N/2)=-1;
    case 2 % DD
      alpha_n=zeros(N,1);
    case 3 % SC
      tao_n=h./mu_n;
      alpha_n=coth(0.5*tao_n)-2./tao_n;
  end % switch FDM
  
  % N rays to trace, each angle has only 1 ray, no ray-spacing
  % n for each angle, and j for FSR region index
  
  phi0_j_old=ones(J,1);
  q_j_n=zeros(J,N);
  for iIterate=1:maxIterate
    for j=1:J
      for n=1:N
        q_j_n(j,n)=(Sig_ss_j(j))*(phi0_j_old(j)-error_ang_j(j))*0.5+Q_MMS_j_n(j,n);
      end
    end
    phi0_j_new=zeros(J,1);


    %% backward direction
    for n=1:N/2
      psi_in=psi_b2_n(n);
      for j=J:-1:1
        temp=abs(mu_n(n))+abs(alpha_n(n))*Sig_t_j(j)*h*0.5;
        psi_out=((temp-0.5*Sig_t_j(j)*h)*psi_in+(q_j_n(j,n))*h)/(temp+0.5*Sig_t_j(j)*h);
        psi_avg=0.5*(1+abs(alpha_n(n)))*psi_out+0.5*(1-abs(alpha_n(n)))*psi_in;
        phi0_j_new(j)=phi0_j_new(j)+psi_avg*weight_n(n);
        psi_in=psi_out;
      end
    end

%% forward direction
    for n=(N/2+1):N
      psi_in=psi_b1_n(n);
      for j=1:J
        temp=mu_n(n)+alpha_n(n)*Sig_t_j(j)*h*0.5;
        psi_out=((temp-0.5*Sig_t_j(j)*h)*psi_in+(q_j_n(j,n))*h)/(temp+0.5*Sig_t_j(j)*h);
        psi_avg=0.5*(1+alpha_n(n))*psi_out+0.5*(1-alpha_n(n))*psi_in;
        phi0_j_new(j)=phi0_j_new(j)+psi_avg*weight_n(n);
        psi_in=psi_out;
      end
    end

%%
    % test for convergence
%     error=norm(phi0_j_new-phi0_j_old);
    error=max(abs(phi0_j_new-phi0_j_old)./(phi0_j_new+delta));
    if error<epsilon_phi0
      break;
    end
    phi0_j_old=phi0_j_new;
  end  

  phi0_j=phi0_j_new;
    
end
