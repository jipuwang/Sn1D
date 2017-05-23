% one region continuous problem
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
function [ ] = Sn_1D_kord_MMSSrc_1Reg()
%     FDMs=[1,2,3];
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
%     nCells=160*[1,2,4,8,16];
%     nCells=[10,20,40,80];
%     nCells=[60];
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
%     Sig_r_g1=zeros(J,1); % not useful.
    Sig_t_g1=zeros(J,1); % total xs.
%     Sig_f_g1=zeros(J,1); % not useful.
    %moderator
    for j=1:J/3
        Sig_ss_g1(j)=0.2020315*multiplier;
        Sig_rm_g1(j)=0.02978*multiplier;
%         Sig_r_g1(j)=0.0002885;
        Sig_t_g1(j)=0.2321;
%         Sig_f_g1(j)=0.0;
    end
    %fuel
    for j=J/3+1:J*2/3
        Sig_ss_g1(j)=0.211762;
        Sig_rm_g1(j)=0.01657;
%         Sig_r_g1(j)=0.009514;
        Sig_t_g1(j)=0.244;
%         Sig_f_g1(j)=0.006154;
    end
    %moderator
    for j=J*2/3+1:J
        Sig_ss_g1(j)=0.2020315*multiplier;
        Sig_rm_g1(j)=0.02978*multiplier;
%         Sig_r_g1(j)=0.0002885;
        Sig_t_g1(j)=0.2321;
%         Sig_f_g1(j)=0.0;
    end

    % group 2, thermal group
    Sig_ss_g2=zeros(J,1); % Self-scattering xs, no up-scat from g2 to g1.
%     Sig_r_g2=zeros(J,1);  
    Sig_t_g2=zeros(J,1);
%     Sig_f_g2=zeros(J,1);
    %moderator
    for j=1:J/3
        Sig_ss_g2(j)=1.090426*multiplier;
%         Sig_r_g2(j)=0.008974;
        Sig_t_g2(j)=1.0994;
%         Sig_f_g2(j)=0.0;
    end
    %fuel
    for j=J/3+1:J*2/3
        Sig_ss_g2(j)=0.6122;
%         Sig_r_g2(j)=0.1026;
        Sig_t_g2(j)=0.8477;
%         Sig_f_g2(j)=0.1329;
    end
    %moderator
    for j=J*2/3+1:J
        Sig_ss_g2(j)=1.090426*multiplier;
%         Sig_r_g2(j)=0.008974;
        Sig_t_g2(j)=1.0994;
%         Sig_f_g2(j)=0.0;
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

%% For the fast group
    % Bounday conditions and source
    % phi=Reg1:y = 1.5593897658139900E-03 x2 -1.3533502228265900E-03 x +6.3255899575709700E-02
    psi_b1_n=0.5*(6.3255899575709700E-02)*ones(N,1); % the -tive half N/2 are not useful
    % phi=Reg3:y = 1.6006807673739400E-03 x2 -1.8963240555487500E-01 x +5.6794357067491300E+00
    Z_Reg1=Z/3;
    psi_b2_n=0.5*(1.5593897658139900E-03*Z_Reg1*Z_Reg1 ...
        -1.3533502228265900E-03*Z_Reg1 ...
        +6.3255899575709700E-02)*ones(N,1); % the +tive half N/2 are not useful

    Q_MMS_n_j_g1=zeros(N,J); % preallocate memory, avg'ed over tau_(j-1/2) and tau_(j+1/2)
    phi0_MMS_j_g1=zeros(J,1);
                                
    % Fast MMS Source in Region 1
    % phi_MMS_1g = ( +1.5593897658139900E-03 x^2 -1.3533502228265900E-03 x +6.3255899575709700E-02 )
    % mu_n/2 * ( +3.1187795316279800E-03 x -1.3533502228265900E-03 ) ...
    % + (Sig_t-Sig_ss)*0.5 * phi_MMS_1g
    % (Sig_t-Sig_ss)*0.5 in this case is 1.5034250000000000E-02
    for j=1:J/3
        phi0_MMS_j_g1(j)= ...
                ( +1.5593897658139900E-03 /3*((j-1)*(j-1)+(j-1)*j+j*j)*h*h ...
                  -1.3533502228265900E-03 *0.5*(j-1+j)*h ...
                  +6.3255899575709700E-02 );
        for n=1:N
        Q_MMS_n_j_g1(n,j)=mu_n(n)*0.5*(+3.1187795316279800E-03 *0.5*(j-1+j)*h -1.3533502228265900E-03) ...
            +(Sig_t_g1(j)-Sig_ss_g1(j))*0.5 * phi0_MMS_j_g1(j);
        end % n
    end % j
    % Fast MMS Source in Region 2
    % phi_MMS_1g = ( -4.4443431764232400E-03 x^2 +2.6590845920514000E-01 x -2.8376294939004300E+00 )
    % mu_n/2 * ( -8.8886863528464800E-03 x +2.6590845920514000E-01 ) ...
    % + (Sig_t-Sig_ss)*0.5 * phi_MMS_1g
    % (Sig_t-Sig_ss)*0.5 in this case is 1.6119000000000000E-02
    for j=J/3+1:2/3*J
        phi0_MMS_j_g1(j)= ...
                ( -4.4443431764232400E-03 /3*((j-1)*(j-1)+(j-1)*j+j*j)*h*h ...
                  +2.6590845920514000E-01 *0.5*(j-1+j)*h ...
                  -2.8376294939004300E+00 );
        for n=1:N
        Q_MMS_n_j_g1(n,j)=mu_n(n)*0.5*(-8.8886863528464800E-03 *0.5*(j-1+j)*h +2.6590845920514000E-01) ...
            +(Sig_t_g1(j)-Sig_ss_g1(j))*0.5 * phi0_MMS_j_g1(j);
        end % n
    end % j
    % Fast MMS Source in Region 3
    % phi_MMS_1g = ( +1.6006807673739400E-03 x^2 -1.8963240555487500E-01 x +5.6794357067491300E+00 )
    % mu_n/2 * ( +3.2013615347478800E-03 x -1.8963240555487500E-01 ) ...
    % + (Sig_t-Sig_ss)*0.5 * phi_MMS_1g
    % (Sig_t-Sig_ss)*0.5 in this case is 1.5034250000000000E-02
    for j=2/3*J+1:J
      phi0_MMS_j_g1(j)= ...
                ( +1.6006807673739400E-03 /3*((j-1)*(j-1)+(j-1)*j+j*j)*h*h ...
                  -1.8963240555487500E-01 *0.5*(j-1+j)*h ...
                  +5.6794357067491300E+00 );
        for n=1:N
        Q_MMS_n_j_g1(n,j)=mu_n(n)*0.5*(+3.2013615347478800E-03 *0.5*(j-1+j)*h -1.8963240555487500E-01) ...
            +(Sig_t_g1(j)-Sig_ss_g1(j))*0.5 * phi0_MMS_j_g1(j);
        end % n
    end % j

%% sacttering source iteration
    S_n_j=zeros(N,J);
    kMax=2000;
    for k=1:kMax
        if k==1
            phi0_g1_j_old=ones(J,1);
        end
        for j=1:J
            S_n_j(:,j)=0.5*(Sig_ss_g1(j)*phi0_g1_j_old(j));
        end % j

        %determine the flux based on this source.
        phi0_j_g1=zeros(J,1);

        for n=(N/2+1):N
            for j=1:J/3
                if j==1
                    psi_n_in=psi_b1_n(n);
                end
                temp=mu_n(n)+alpha_n(n)*Sig_t_g1(j)*h*0.5;
                psi_n_out=((temp-0.5*Sig_t_g1(j)*h)*psi_n_in+(S_n_j(n,j)+Q_MMS_n_j_g1(n,j))*h)/(temp+0.5*Sig_t_g1(j)*h);
                psi_n_avg=0.5*(1+alpha_n(n))*psi_n_out+0.5*(1-alpha_n(n))*psi_n_in;
                phi0_j_g1(j)=phi0_j_g1(j)+psi_n_avg*weight_n(n);
                psi_n_in=psi_n_out;
            end
        end

        for n=1:N/2
            for j=J/3:-1:1
                if j==J/3
                    psi_n_in=psi_b2_n(n);
                end
                temp=abs(mu_n(n))+abs(alpha_n(n))*Sig_t_g1(j)*h*0.5;
                psi_n_out=((temp-0.5*Sig_t_g1(j)*h)*psi_n_in+(S_n_j(n,j)+Q_MMS_n_j_g1(n,j))*h)/(temp+0.5*Sig_t_g1(j)*h);
                psi_n_avg=0.5*(1+abs(alpha_n(n)))*psi_n_out+0.5*(1-abs(alpha_n(n)))*psi_n_in;
                phi0_j_g1(j)=phi0_j_g1(j)+psi_n_avg*weight_n(n);
                psi_n_in=psi_n_out;
            end
        end

        error=max(abs(phi0_j_g1(1:J/3)-phi0_g1_j_old(1:J/3))./(phi0_j_g1(1:J/3)+delta))
        if error<epsilon
            break;
        else
            phi0_g1_j_old(1:J/3)=phi0_j_g1(1:J/3);
        end
    end
    k
    if k>kMax-1
        display 'Not converging';
    end
    figure(17);plot(phi0_j_g1,'-*');
    title('fast flux');

%% For thermal group
    % Bounday conditions and source
    % phi=Reg1:y = +3.7420664087460500E-04 x2 +3.2444814149334100E-03 x +1.2139249518777500E-02
    psi_b1_n=0.5*1.2139249518777500E-02*ones(N,1); % the -tive half N/2 are not useful
    % phi=Reg3:y = +3.6321032591887700E-04x2 -4.7314818437341200E-02x +1.5394028032516800E+00
    psi_b2_n=0.5*(+3.7420664087460500E-04*Z_Reg1*Z_Reg1 ...
                  +3.2444814149334100E-03*Z_Reg1 ...
                  +1.2139249518777500E-02)*ones(N,1); % the +tive half N/2 are not useful
    Q_n_j_g2=zeros(N,J); % Inscatter/downscatter from group 1. 
    for j=1:J
        Q_n_j_g2(:,j)=0.5*Sig_rm_g1(j)*phi0_j_g1(j);
    end % j

    Q_MMS_n_j_g2=zeros(N,J);
    phi0_MMS_j_g2=zeros(J,1);
    % Thermal MMS Source in Region 1
    % phi_MMS_2g = ( +3.7420664087460500E-04 x^2 +3.2444814149334100E-03 x +1.2139249518777500E-02 )
    % mu_n/2 * ( +7.4841328174921000E-04 x +3.2444814149334100E-03 ) ...
    % + (Sig_t-Sig_ss)*0.5 * phi_MMS_g2
    % - Sig_g1_rm*0.5 * (phi_g1_MMS)
    % (Sig_t-Sig_ss)*0.5* in this case is +4.4870000000000200E-03
    for j=1:J/3
        phi0_MMS_j_g2(j)= ...
                ( +3.7420664087460500E-04 /3*((j-1)*(j-1)+(j-1)*j+j*j)*h*h ...
                  +3.2444814149334100E-03 *0.5*(j-1+j)*h ...
                  +1.2139249518777500E-02 );
        for n=1:N
        Q_MMS_n_j_g2(n,j)=mu_n(n)*0.5*(+7.4841328174921000E-04 *0.5*(j-1+j)*h +3.2444814149334100E-03) ...
            +(Sig_t_g2(j)-Sig_ss_g2(j))*0.5 * phi0_MMS_j_g2(j) ...
            -Sig_rm_g1(j)*0.5 * phi0_MMS_j_g1(j);
        end % n
    end % j
    % Thermal MMS Source in Region 2
    % phi_MMS_g2 = ( -1.2351628643713800E-03 x^2 +7.3932137536755800E-02 x -7.5807942857102300E-01 )
    % mu_n/2 * ( -2.4703257287427600E-03 x +7.3932137536755800E-02 ) ...
    % + (Sig_t-Sig_ss)*0.5 * phi_MMS_g2
    % - Sig_g1_rm*0.5 * (phi_MMS_g1)
    % (Sig_t-Sig_ss)*0.5 in this case is +1.1775000000000000E-01
    for j=J/3+1:2/3*J
        phi0_MMS_j_g2(j)= ...
                ( -1.2351628643713800E-03 /3*((j-1)*(j-1)+(j-1)*j+j*j)*h*h ...
                  +7.3932137536755800E-02 *0.5*(j-1+j)*h ...
                  -7.5807942857102300E-01 );
        for n=1:N
        Q_MMS_n_j_g2(n,j)=mu_n(n)*0.5*(-2.4703257287427600E-03 *0.5*(j-1+j)*h +7.3932137536755800E-02) ...
            +(Sig_t_g2(j)-Sig_ss_g2(j))*0.5 * phi0_MMS_j_g2(j) ...
            -Sig_rm_g1(j)*0.5 * phi0_MMS_j_g1(j);
        end % n
    end % j

    % Thermal MMS Source in Region 3
    % phi_MMS_g2 = ( +3.6321032591887700E-04 x^2 -4.7314818437341200E-02 x +1.5394028032516800E+00 )
    % mu_n/2* ( +7.2642065183775400E-04 x -4.7314818437341200E-02 ) ...
    % + (Sig_t-Sig_ss)*0.5* phi_MMS_g2
    % - Sig_g1_rm*0.5 * (phi_MMS_g1)
    % (Sig_t-Sig_ss)*0.5 in this case is +4.4870000000000200E-03
    for j=2/3*J+1:J
        phi0_MMS_j_g2(j)= ...
                ( +3.6321032591887700E-04 /3*((j-1)*(j-1)+(j-1)*j+j*j)*h*h ...
                  -4.7314818437341200E-02 *0.5*(j-1+j)*h ...
                  +1.5394028032516800E+00 );
        for n=1:N
        Q_MMS_n_j_g2(n,j)=mu_n(n)*0.5*(+7.2642065183775400E-04 *0.5*(j-1+j)*h -4.7314818437341200E-02) ...
            +(Sig_t_g2(j)-Sig_ss_g2(j))*0.5 * phi0_MMS_j_g2(j) ...
            -Sig_rm_g1(j)*0.5 * phi0_MMS_j_g1(j);
        end % n
    end % j

    for k=1:kMax
        if k==1
            phi0_j_g2_old=ones(J,1);
        end
        for j=1:J
            S_n_j(:,j)=0.5*(Sig_ss_g2(j)*phi0_j_g2_old(j));
        end % j
        %determine the flux based on this source.
        phi0_j_g2=zeros(J,1);

        for n=(N/2+1):N
            for j=1:J/3
                if j==1
                    psi_n_in=psi_b1_n(n);
                end
                temp=mu_n(n)+alpha_n(n)*Sig_t_g2(j)*h*0.5;
                psi_n_out=((temp-0.5*Sig_t_g2(j)*h)*psi_n_in+(S_n_j(n,j)+Q_n_j_g2(n,j)+Q_MMS_n_j_g2(n,j))*h)/(temp+0.5*Sig_t_g2(j)*h);
                psi_n_avg=0.5*(1+alpha_n(n))*psi_n_out+0.5*(1-alpha_n(n))*psi_n_in;
                phi0_j_g2(j)=phi0_j_g2(j)+psi_n_avg*weight_n(n);
                psi_n_in=psi_n_out;
            end
        end

        for n=1:N/2
            for j=J/3:-1:1
                if j==J/3
                    psi_n_in=psi_b2_n(n);
                end
                temp=abs(mu_n(n))+abs(alpha_n(n))*Sig_t_g2(j)*h*0.5;
                psi_n_out=((temp-0.5*Sig_t_g2(j)*h)*psi_n_in+(S_n_j(n,j)+Q_n_j_g2(n,j)+Q_MMS_n_j_g2(n,j))*h)/(temp+0.5*Sig_t_g2(j)*h);
                psi_n_avg=0.5*(1+abs(alpha_n(n)))*psi_n_out+0.5*(1-abs(alpha_n(n)))*psi_n_in;
                phi0_j_g2(j)=phi0_j_g2(j)+psi_n_avg*weight_n(n);
                psi_n_in=psi_n_out;
            end
        end

        error=max(abs(phi0_j_g2(1:J/3)-phi0_j_g2_old(1:J/3))./(phi0_j_g2(1:J/3)+delta))
        if error<epsilon
            break;
        else
            phi0_j_g2_old(1:J/3)=phi0_j_g2(1:J/3);
        end
    end
    k
    if k>kMax-1
        display 'Not converging';
    end
    % plots
    figure(18);plot(phi0_j_g2,'-*'); title('thermal flux');
    figure(19);plot(phi0_j_g1,'-*'); hold on; plot(phi0_j_g2,'-*');
    title('with Manufactured Source');
    legend('fast flux','thermal flux'); hold off;
    
    % getter for regression test
    global getter_phi0_j_g1 getter_phi0_j_g2
    getter_phi0_j_g1=phi0_j_g1;
    getter_phi0_j_g2=phi0_j_g2;
    
    % grid function norm
    error_j_g1=norm(phi0_j_g1(1:J/3)-phi0_MMS_j_g1(1:J/3))/sqrt(J);
    error_j_g2=norm(phi0_j_g2(1:J/3)-phi0_MMS_j_g2(1:J/3))/sqrt(J);
end

