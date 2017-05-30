%% Sn 1D solver - with fixed source
% This solver solves planar-geometry Sn equation.
% Spatial: 1D, 3 different FDMs variable z.
% Step Differencing Method (Step), Diamond Differencing (DD), and Step
% Characteristics (SC).
% Multiple region  problem are defined by the cross section.
% E: 2 energy groups.
% Angle: Sn.
% Time: steady state.
% MMS Source: Fixed specified source, angle-dependent
%% Main Function
function [ ] = heavyBoron_MMSSrc_cont_flux()
%   FDMs=[1,2,3];
  FDMs=[2];
  % 1:Step; 2:DD; 3:SC;
  nFDMs=size(FDMs,2);
  for iFDM=1:nFDMs
    FDM=FDMs(iFDM);
    Sn_1D_spatial(FDM);
  end
end
%% Grid Refinement Study
function [ ] = Sn_1D_spatial(FDM)
%   nCells=160*[1,2,4,8,16];
%   nCells=[10,20,40,80];
%   nCells=[60];
  nCells=[1,2,4,8]*60*2%*4*4*4;
  nRuns=size(nCells,2);
  RMSs_phi0_g1=zeros(nRuns,1);
  RMSs_phi0_g2=zeros(nRuns,1);
  for iRun=1:nRuns
    nCell=nCells(iRun);
    [RMSs_phi0_g1(iRun),RMSs_phi0_g2(iRun)]=Sn_1D_RMS(nCell,FDM);
    % process error_j_g1,error_j_g2
  end %iRun

  % display the results
  if iRun==nRuns
    switch FDM
      case 1
        disp('FDM: Step Differencing');
      case 2
        disp('FDM: Diamond Difference');
      case 3
        disp('FDM: Step Characteristics');
    end
    display(RMSs_phi0_g1);
    display(RMSs_phi0_g2);
  end %if

  p_phi0_g1=zeros(nRuns-1,1);
  p_phi0_g2=zeros(nRuns-1,1);
  for i=1:nRuns-1
    p_phi0_g1(i)=log(RMSs_phi0_g1(i)/RMSs_phi0_g1(i+1))/log(nCells(i+1)/nCells(i));
    p_phi0_g2(i)=log(RMSs_phi0_g2(i)/RMSs_phi0_g2(i+1))/log(nCells(i+1)/nCells(i));
  end

  display(p_phi0_g1);
  display(p_phi0_g2);

end

%% Error Estimate - RMS
function [error_j_g1,error_j_g2]=Sn_1D_RMS(J,FDM)
  % J: num. of spatial cells
  % Z: system width [cm]
  Z=60;
  % cell width [cm]
  h=Z/J;
  % Quadrature Set
  global N; N=4;%16;
  [mu_n,weight_n]=lgwt(N,-1,1);mu_n=flipud(mu_n);
  % convergence criterion
  epsilon=1E-10;
  delta=1E-13;

  % specify material cross sections.
  % scattering xs multiplier 0<multiplier<1.008
  multiplier=1*1.0;
  % group 1, fast group
  Sig_ss_g1=zeros(J,1); % Self-scattering xs.
  Sig_rm_g1=zeros(J,1); % Removal xs from g1 to g2.
%   Sig_r_g1=zeros(J,1); % not useful.
  Sig_t_g1=zeros(J,1); % total xs.
%   Sig_f_g1=zeros(J,1); % not useful.
  %moderator
  for j=1:J/3
    Sig_ss_g1(j)=0.2020315*multiplier;
    Sig_rm_g1(j)=0.02978*multiplier;
%     Sig_r_g1(j)=0.0002885;
    Sig_t_g1(j)=0.2321;
%     Sig_f_g1(j)=0.0;
  end
  %fuel
  for j=J/3+1:J*2/3
    Sig_ss_g1(j)=0.211762;
    Sig_rm_g1(j)=0.01657;
%     Sig_r_g1(j)=0.009514;
    Sig_t_g1(j)=0.244;
%     Sig_f_g1(j)=0.006154;
  end
  %moderator
  for j=J*2/3+1:J
    Sig_ss_g1(j)=0.2020315*multiplier;
    Sig_rm_g1(j)=0.02978*multiplier;
%     Sig_r_g1(j)=0.0002885;
    Sig_t_g1(j)=0.2321;
%     Sig_f_g1(j)=0.0;
  end

  % group 2, thermal group
  Sig_ss_g2=zeros(J,1); % Self-scattering xs, no up-scat from g2 to g1.
%   Sig_r_g2=zeros(J,1);
  Sig_t_g2=zeros(J,1);
%   Sig_f_g2=zeros(J,1);
  %moderator
  for j=1:J/3
    Sig_ss_g2(j)=1.090426*multiplier;
%     Sig_r_g2(j)=0.008974;
    Sig_t_g2(j)=1.0994;
%     Sig_f_g2(j)=0.0;
  end
  %fuel
  for j=J/3+1:J*2/3
    Sig_ss_g2(j)=0.6122;
%     Sig_r_g2(j)=0.1026;
    Sig_t_g2(j)=0.8477;
%     Sig_f_g2(j)=0.1329;
  end
  %moderator
  for j=J*2/3+1:J
    Sig_ss_g2(j)=1.090426*multiplier;
%     Sig_r_g2(j)=0.008974;
    Sig_t_g2(j)=1.0994;
%     Sig_f_g2(j)=0.0;
  end

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
  % Bounday conditions and source
  % reg1 expression evaluated at x=0 and take half
  psi_b1_n=0.5*reg1(0)*ones(N,1); % the -tive half N/2 are not useful
  % reg3 expression evaluated at x=Z and take half
  psi_b2_n=0.5*reg3(Z)*ones(N,1); % the +tive half N/2 are not useful

  Q_MMS_n_j_g1=zeros(N,J); % preallocate memory, avg'ed over tau_(j-1/2) and tau_(j+1/2)
  % MMS source: mu_n/2 * derivative(phi_MMS_g1) ...
  % + (Sig_t-Sig_ss)*0.5 * phi_MMS_g1
  phi0_MMS_j_g1=zeros(J,1);

  % Fast MMS Source in Region 1
  % phi_MMS_g1 = reg1
  for j=1:J/3
    x_L=(j-1)*h;x_R=j*h;
    phi0_MMS_j_g1(j)= 1/h*integral(reg1,x_L,x_R);
    for n=1:N
    Q_MMS_n_j_g1(n,j)=mu_n(n)*0.5*(2*reg1x2 *0.5*(j-1+j)*h +reg1x1) ...
      +(Sig_t_g1(j)-Sig_ss_g1(j))*0.5 * phi0_MMS_j_g1(j);
    end % n
  end % j
  % Fast MMS Source in Region 2
  % For the next two regions, nothing changed exception for the coeff.
  % phi_MMS_g1 = reg2
  for j=J/3+1:2/3*J
    x_L=(j-1)*h;x_R=j*h;
    phi0_MMS_j_g1(j)= 1/h*integral(reg2,x_L,x_R);
    for n=1:N
    Q_MMS_n_j_g1(n,j)=mu_n(n)*0.5*(2*reg2x2 *0.5*(j-1+j)*h +reg2x1) ...
      +(Sig_t_g1(j)-Sig_ss_g1(j))*0.5 * phi0_MMS_j_g1(j);
    end % n
  end % j
  % Fast MMS Source in Region 3
  % phi_MMS_g1 = reg3
  for j=2/3*J+1:J
    x_L=(j-1)*h;x_R=j*h;
    phi0_MMS_j_g1(j)= 1/h*integral(reg3,x_L,x_R);
    for n=1:N
    Q_MMS_n_j_g1(n,j)=mu_n(n)*0.5*(2*reg3x2 *0.5*(j-1+j)*h +reg3x1) ...
      +(Sig_t_g1(j)-Sig_ss_g1(j))*0.5 * phi0_MMS_j_g1(j);
    end % n
  end % j

%% Non-MMS stuff
  % Start converging self-scattering
  self_scat_n_j=zeros(N,J);
  kMax=2000;
  for k=1:kMax
    if k==1
      phi0_g1_j_old=ones(J,1);
    end
    for j=1:J
      self_scat_n_j(:,j)=0.5*(Sig_ss_g1(j)*phi0_g1_j_old(j));
    end % j

    %determine the flux based on this source.
    phi0_j_g1=zeros(J,1);

    for n=(N/2+1):N
      for j=1:J
        if j==1
          psi_n_in=psi_b1_n(n);
        end
        temp=mu_n(n)+alpha_n(n)*Sig_t_g1(j)*h*0.5;
        psi_n_out=((temp-0.5*Sig_t_g1(j)*h)*psi_n_in+(self_scat_n_j(n,j)+Q_MMS_n_j_g1(n,j))*h)/(temp+0.5*Sig_t_g1(j)*h);
        psi_n_avg=0.5*(1+alpha_n(n))*psi_n_out+0.5*(1-alpha_n(n))*psi_n_in;
        phi0_j_g1(j)=phi0_j_g1(j)+psi_n_avg*weight_n(n);
        psi_n_in=psi_n_out;
      end
    end

    for n=1:N/2
      for j=J:-1:1
        if j==J
          psi_n_in=psi_b2_n(n);
        end
        temp=abs(mu_n(n))+abs(alpha_n(n))*Sig_t_g1(j)*h*0.5;
        psi_n_out=((temp-0.5*Sig_t_g1(j)*h)*psi_n_in+(self_scat_n_j(n,j)+Q_MMS_n_j_g1(n,j))*h)/(temp+0.5*Sig_t_g1(j)*h);
        psi_n_avg=0.5*(1+abs(alpha_n(n)))*psi_n_out+0.5*(1-abs(alpha_n(n)))*psi_n_in;
        phi0_j_g1(j)=phi0_j_g1(j)+psi_n_avg*weight_n(n);
        psi_n_in=psi_n_out;
      end
    end

    error=max(abs(phi0_j_g1-phi0_g1_j_old)./(phi0_j_g1+delta))
    if error<epsilon
      break;
    else
      phi0_g1_j_old=phi0_j_g1;
    end
  end
  k
  if k>kMax-1
    display 'Not converging';
  end
  figure(13);plot(phi0_j_g1,'-*');
  title('fast flux');

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
  % Bounday conditions and source
  % reg1 expression evaluated at x=0 and take half
  psi_b1_n=0.5*reg1(0)*ones(N,1); % the -tive half N/2 are not useful
  % reg3 expression evaluated at x=Z and take half
  psi_b2_n=0.5*reg3(Z)*ones(N,1); % the +tive half N/2 are not useful


  Q_MMS_n_j_g2=zeros(N,J);
  % MMS source: mu_n/2 * derivative(phi_MMS_g2) ...
  % + (Sig_t-Sig_ss)*0.5 * phi_MMS_g2
  % - Sig_g1_rm*0.5 * (phi_MMS_g1)
  phi0_MMS_j_g2=zeros(J,1);
  % Thermal MMS Source in Region 1
  % phi_MMS_g1 = reg1
  for j=1:J/3
    x_L=(j-1)*h;x_R=j*h;
    phi0_MMS_j_g2(j)= 1/h*integral(reg1,x_L,x_R);
    for n=1:N
    Q_MMS_n_j_g2(n,j)=mu_n(n)*0.5*(2*reg1x2 *0.5*(j-1+j)*h +reg1x1) ...
      +(Sig_t_g2(j)-Sig_ss_g2(j))*0.5 * phi0_MMS_j_g2(j) ...
      -Sig_rm_g1(j)*0.5 * phi0_MMS_j_g1(j);
    end % n
  end % j
  % Thermal MMS Source in Region 2
  % phi_MMS_g1 = reg2
  % The only difference in the following 2 regions are
  % reg2x2,reg2x1,reg2x0
  for j=J/3+1:2/3*J
    x_L=(j-1)*h;x_R=j*h;
    phi0_MMS_j_g2(j)= 1/h*integral(reg2,x_L,x_R);
    for n=1:N
    Q_MMS_n_j_g2(n,j)=mu_n(n)*0.5*(2*reg2x2 *0.5*(j-1+j)*h +reg2x1) ...
      +(Sig_t_g2(j)-Sig_ss_g2(j))*0.5 * phi0_MMS_j_g2(j) ...
      -Sig_rm_g1(j)*0.5 * phi0_MMS_j_g1(j);
    end % n
  end % j

  % Thermal MMS Source in Region 3
  % phi_MMS_g1 = reg3
  for j=2/3*J+1:J
    x_L=(j-1)*h;x_R=j*h;
    phi0_MMS_j_g2(j)= 1/h*integral(reg3,x_L,x_R);
    for n=1:N
    Q_MMS_n_j_g2(n,j)=mu_n(n)*0.5*(2*reg3x2 *0.5*(j-1+j)*h +reg3x1) ...
      +(Sig_t_g2(j)-Sig_ss_g2(j))*0.5 * phi0_MMS_j_g2(j) ...
      -Sig_rm_g1(j)*0.5 * phi0_MMS_j_g1(j);
    end % n
  end % j

%% Non-MMS stuff
  % Inscattering source component
  inscat_n_j_g2=zeros(N,J); % Inscatter/downscatter from group 1.
  for j=1:J
    inscat_n_j_g2(:,j)=0.5*Sig_rm_g1(j)*phi0_j_g1(j);
  end % j

  % Start converging self-scattering
  self_scat_n_j=zeros(N,J);
  for k=1:kMax
    if k==1
      phi0_j_g2_old=ones(J,1);
    end
    for j=1:J
      self_scat_n_j(:,j)=0.5*(Sig_ss_g2(j)*phi0_j_g2_old(j));
    end % j
    %determine the flux based on this source.
    phi0_j_g2=zeros(J,1);

    for n=(N/2+1):N
      for j=1:J
        if j==1
          psi_n_in=psi_b1_n(n);
        end
        temp=mu_n(n)+alpha_n(n)*Sig_t_g2(j)*h*0.5;
        psi_n_out=((temp-0.5*Sig_t_g2(j)*h)*psi_n_in+(self_scat_n_j(n,j)+inscat_n_j_g2(n,j)+Q_MMS_n_j_g2(n,j))*h)/(temp+0.5*Sig_t_g2(j)*h);
        psi_n_avg=0.5*(1+alpha_n(n))*psi_n_out+0.5*(1-alpha_n(n))*psi_n_in;
        phi0_j_g2(j)=phi0_j_g2(j)+psi_n_avg*weight_n(n);
        psi_n_in=psi_n_out;
      end
    end

    for n=1:N/2
      for j=J:-1:1
        if j==J
          psi_n_in=psi_b2_n(n);
        end
        temp=abs(mu_n(n))+abs(alpha_n(n))*Sig_t_g2(j)*h*0.5;
        psi_n_out=((temp-0.5*Sig_t_g2(j)*h)*psi_n_in+(self_scat_n_j(n,j)+inscat_n_j_g2(n,j)+Q_MMS_n_j_g2(n,j))*h)/(temp+0.5*Sig_t_g2(j)*h);
        psi_n_avg=0.5*(1+abs(alpha_n(n)))*psi_n_out+0.5*(1-abs(alpha_n(n)))*psi_n_in;
        phi0_j_g2(j)=phi0_j_g2(j)+psi_n_avg*weight_n(n);
        psi_n_in=psi_n_out;
      end
    end

    error=max(abs(phi0_j_g2-phi0_j_g2_old)./(phi0_j_g2+delta))
    if error<epsilon
      break;
    else
      phi0_j_g2_old=phi0_j_g2;
    end
  end
  k
  if k>kMax-1
    display 'Not converging';
  end
  % plots
  figure(14);plot(phi0_j_g2,'-*'); title('thermal flux');
  figure(15);plot(phi0_j_g1,'-*'); hold on; plot(phi0_j_g2,'-*');
  title('with Continuous Manufactured Source');
  legend('fast flux','thermal flux'); hold off;

  % getter for regression test
  global getter_phi0_j_g1 getter_phi0_j_g2
  getter_phi0_j_g1=phi0_j_g1;
  getter_phi0_j_g2=phi0_j_g2;

  % grid function norm
  error_j_g1=norm(phi0_j_g1-phi0_MMS_j_g1)/sqrt(J);
  error_j_g2=norm(phi0_j_g2-phi0_MMS_j_g2)/sqrt(J);
end

