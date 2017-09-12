% A sine and sine solution generator including
    % Discretized analytical solution
    % Manufactured boundary conditions
    % Manufactured source
function [phi0_MMS_j,psi_b1_n,psi_b2_n,Q_MMS_j_n,error_ang_j,phi0_guess_j,k_guess]=...
          manufacturer_SnEig(J,N,Tau,mat,assumedSoln,k_MMS)
  % input parameters
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
  if ~exist('assumedSoln','var')
    assumedSoln='constant';
  end
  if ~exist('k_MMS','var')
    k_MMS=1.02;
  end
  
  % Material
  Sig_t_j=mat.Sig_t_j;
  Sig_ss_j=mat.Sig_ss_j;
  Sig_gamma_j=mat.Sig_gamma_j;
  Sig_f_j=mat.Sig_f_j;
  nuSig_f_j=mat.nuSig_f_j;
  kappaSig_f_j=mat.kappaSig_f_j;
  k_F=mat.thermal_cond_k_j(1);

  h=Tau/J;
  [mu_n,weight_n]=lgwt(N,-1,1); mu_n=flipud(mu_n);
  
  %% Manufactured Solutions for both fields
  % Options includes: assumedSoln='constant'; 'linear'; 'quadratic';
  % 'plus1Sqrt'; 'other_anisotropic';

  switch(assumedSoln)
    case('constant')
      % Manufactured neutronics solution \psi(x,\mu)=sin(pi*x/Tau), 0<x<Tau
      psi_MMS =@(x,mu) (1.0+0.0*x).*(1.0+0.0*mu);
      psi_MMS_Diff =@(x,mu) (0.0+0.0*x).*(1.0+0.0*mu);
    case('linear')
      psi_MMS =@(x,mu) 1.0+x.*exp(mu);
      psi_MMS_Diff =@(x,mu) (1.0+x*0.0).*exp(mu);
    case('quadratic')
      % Manufactured neutronics solution \psi(x,\mu)=1.0, 0<x<Tau
      psi_MMS =@(x,mu) 1.0+x.*x.*exp(mu);
      psi_MMS_Diff =@(x,mu) (2*x).*exp(mu);
    case('plus1Sqrt')
      % Manufactured neutronics solution \psi(x,\mu)=1.0, 0<x<Tau
      psi_MMS =@(x,mu) sqrt(x+1).*(1.0+0.0*mu);
      psi_MMS_Diff =@(x,mu) 0.5./sqrt(x+1).*(1.0+0.0*mu);
    case('flat_expMu')
      psi_MMS =@(x,mu) (1.0+0.0*x).*exp(mu);
      psi_MMS_Diff =@(x,mu) (0.0+0.0*x).*exp(mu);
    case('sine')
      psi_MMS =@(x,mu) sin(pi*x/Tau).*exp(mu);
      psi_MMS_Diff =@(x,mu) pi/Tau*cos(pi*x/Tau).*exp(mu);
  end
  
  Sig_gamma =@(x) Sig_gamma_j(1)+0.0*x;
  Sig_ss =@(x) Sig_ss_j(1)+x*0;
  Sig_f =@(x) Sig_f_j(1)+x*0;
  nuSig_f =@(x) nuSig_f_j(1)+x*0;
  Sig_t =@(x) Sig_ss(x)+Sig_f(x)+Sig_gamma(x);
  
%% Manufactured scalar flux and source
  phi0_MMS =@(x) integral(@(mu) psi_MMS(x,mu), -1,1);
  % MMS source: mu_n * derivative(psi_MMS) +Sig_t* psi_MMS ...
  % -(Sig_ss+nuSig_f/k)*0.5*phi0_MMS;
  Q_MMS =@(x,mu) mu*psi_MMS_Diff(x,mu) +Sig_t(x).*psi_MMS(x,mu) ...
    -(Sig_ss(x)+nuSig_f(x)/k_MMS)*0.5.*phi0_MMS(x);
  
  %% For MoC MMS solution and problem
  % Boundary condition and source
  psi_b1_n=zeros(N,1);
  % psi expression evaluated at x=0
  for n=1:N
      psi_b1_n(n)=psi_MMS(0,mu_n(n)); % n=N/2+1:N % mu>0
  end
  % psi expression evaluated at x=Tau
  psi_b2_n=zeros(N,1);
  for n=1:N
      psi_b2_n(n)=psi_MMS(Tau,mu_n(n)); % n=N/2+1:N % mu>0
  end
  
  phi0_MMS_j=zeros(J,1);
  Q_MMS_j_n=zeros(J,N);
  error_ang_j=ones(J,1);

  for j=1:J
    x_L=(j-1)*h;x_R=j*h;
    phi0_MMS_j(j)=1/h*integral(phi0_MMS,x_L,x_R, 'ArrayValued',true);
    numSum=0;
    for n=1:N
        Q_MMS_j_n(j,n)= ...
          1/h*integral(@(x) Q_MMS(x,mu_n(n)),x_L,x_R, 'ArrayValued',true);
        spatialAvg=1/h*integral(@(x) psi_MMS(x,mu_n(n)),x_L,x_R);
        numSum=numSum+weight_n(n)*spatialAvg;
    end % n
    error_ang_j(j)=numSum-phi0_MMS_j(j);
  end % j

  % Determine the phi0_guess_j and k_guess
  phi0_guess_j=ones(J,1);
  k_guess=k_MMS*(sum(nuSig_f_j.*phi0_guess_j)*h)/(sum(nuSig_f_j.*phi0_MMS_j)*h);
end
