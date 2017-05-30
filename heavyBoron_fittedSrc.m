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
function [ ] = Sn_1D_kord_fittedSrc()
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
    nCells=[60]*2%*4*4*4;
    nRuns=size(nCells,2);
    for iRun=1:nRuns
        nCell=nCells(iRun);
        [phi0_g1_j,phi0_g2_j]=Sn_1D_RMS(nCell,FDM);
        % process phi0_g1_j and phi0_g2_j
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
%         display(phi0_g1_j);
%         display(phi0_g2_j);
    end %if
end

%% Error Estimate - RMS
function [phi0_j_g1,phi0_j_g2]=Sn_1D_RMS(J,FDM)
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
% Bounday conditions and source for fast group
    psi_b1_n_g1=zeros(N,1);  % the -tive half N/2 are not useful
    psi_b2_n_g1=zeros(N,1);  % the +tive half N/2 are not useful
%% fitted source
    Q_n_j_g1=zeros(N,J); % fast source
    % Only fuel section has fast fixed source
    % The following fitted source is the fission source built from the
    % flux digitized from the plots in Kord's slides and then used fission 
    % xs to build the formula. 
    sourceFitting = @(x) -5.3033745077014600E-06*x.^5 ...
    +6.8115366378200500E-04*x.^4 -3.4950409538850600E-02*x.^3 ...
    +8.9576122737953600E-01*x.^2 -1.1468981660532800E+01*x ...
    +5.8829465224422900E+01;
    fittedSource = zeros(J/3,1);
    for j=1:J/3/2
        x_L=Z/3+(j-1)*h;
        x_R=x_L+h;
        temp=1/h*(integral(sourceFitting, x_L, x_R));
        fittedSource(j)=temp;
        fittedSource(J/3-j+1)=temp;
    end
    for j=J/3+1:2/3*J % Only the fuel region has this source.
%         Q_nj_1(:,j)=0.5*0.16;
        Q_n_j_g1(:,j)=0.5*fittedSource(j-(J/3));
    end % j
% %%  temporary boundary definition and source. For constant source case.
%     psi_b1_n_g1=zeros(N,1); % the -tive half N/2 are not useful
%     psi_b2_n_g1=zeros(N,1); % the +tive half N/2 are not useful
% 
%     Q_n_j_g1=zeros(N,J);
%     for j=J/3+1:2/3*J % Only the fuel region has fast source. Overwritten.
%         Q_n_j_g1(:,j)=0.5*0.16;
%     end % j
%% This part is for converging the scattering source
    S_n_j=zeros(N,J);
    iterMax=2000;
    for iter=1:iterMax
        if iter==1
            phi0_j_g1_old=ones(J,1);
        end
        % update the scattering source
        for j=1:J
            S_n_j(:,j)=0.5*(Sig_ss_g1(j)*phi0_j_g1_old(j));
        end % j

        phi0_j_g1=zeros(J,1);
        for n=(N/2+1):N
            for j=1:J
                if j==1
                    psi_n_in=psi_b1_n_g1(n);
                end
                temp=mu_n(n)+alpha_n(n)*Sig_t_g1(j)*h*0.5;
                psi_n_out=((temp-0.5*Sig_t_g1(j)*h)*psi_n_in+(S_n_j(n,j)+Q_n_j_g1(n,j))*h)/(temp+0.5*Sig_t_g1(j)*h);
                psi_n_avg=0.5*(1+alpha_n(n))*psi_n_out+0.5*(1-alpha_n(n))*psi_n_in;
                phi0_j_g1(j)=phi0_j_g1(j)+psi_n_avg*weight_n(n);
                psi_n_in=psi_n_out;
                % for debugging
%                 figure(11);plot(phi0_j_g1,'-*');
            end
        end

        for n=1:N/2
            for j=J:-1:1
                if j==J
                    psi_n_in=psi_b2_n_g1(n);
                end
                temp=abs(mu_n(n))+abs(alpha_n(n))*Sig_t_g1(j)*h*0.5;
                psi_n_out=((temp-0.5*Sig_t_g1(j)*h)*psi_n_in+(S_n_j(n,j)+Q_n_j_g1(n,j))*h)/(temp+0.5*Sig_t_g1(j)*h);
                psi_n_avg=0.5*(1+abs(alpha_n(n)))*psi_n_out+0.5*(1-abs(alpha_n(n)))*psi_n_in;
                phi0_j_g1(j)=phi0_j_g1(j)+psi_n_avg*weight_n(n);
                psi_n_in=psi_n_out;
                % for debugging
%                 figure(12);plot(phi0_j_g1,'-*');
            end
        end

        error=max(abs(phi0_j_g1-phi0_j_g1_old)./(phi0_j_g1+delta))
        if error<epsilon
            break;
        else
            phi0_j_g1_old=phi0_j_g1;
        end
    end % iter
    iter
    if iter>iterMax-1
        disp 'Not converging';
    end
    figure(11);plot(phi0_j_g1,'-*');
    title('fast flux');

%% For the thermal group
% Bounday conditions and source for thermal group
    psi_b1_n_2g=zeros(N,1); % the -tive half N/2 are not useful
    psi_b2_n_2g=zeros(N,1); % the +tive half N/2 are not useful
    Q_n_j_2g=zeros(N,J); % no fixed thermal source other than in-scat.
    for j=1:J
        Q_n_j_2g(:,j)=0.5*Sig_rm_g1(j)*phi0_j_g1(j);
    end % j

    for iter=1:iterMax
        if iter==1
            phi0_j_g2_old=ones(J,1);
        end
        % update self-scat source
        for j=1:J
            S_n_j(:,j)=0.5*(Sig_ss_g2(j)*phi0_j_g2_old(j));
        end % j
        phi0_j_g2=zeros(J,1);

        for n=(N/2+1):N
            for j=1:J
                if j==1
                    psi_n_in=psi_b1_n_2g(n);
                end
                temp=mu_n(n)+alpha_n(n)*Sig_t_g2(j)*h*0.5;
                psi_n_out=((temp-0.5*Sig_t_g2(j)*h)*psi_n_in+(S_n_j(n,j)+Q_n_j_2g(n,j))*h)/(temp+0.5*Sig_t_g2(j)*h);
                psi_n_avg=0.5*(1+alpha_n(n))*psi_n_out+0.5*(1-alpha_n(n))*psi_n_in;
                phi0_j_g2(j)=phi0_j_g2(j)+psi_n_avg*weight_n(n);
                psi_n_in=psi_n_out;
            end
        end

        for n=1:N/2
            for j=J:-1:1
                if j==J
                    psi_n_in=psi_b2_n_2g(n);
                end
                temp=abs(mu_n(n))+abs(alpha_n(n))*Sig_t_g2(j)*h*0.5;
                psi_n_out=((temp-0.5*Sig_t_g2(j)*h)*psi_n_in+(S_n_j(n,j)+Q_n_j_2g(n,j))*h)/(temp+0.5*Sig_t_g2(j)*h);
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
    iter
    if iter>iterMax-1
        display 'Not converging';
    end
    
    figure(12);plot(phi0_j_g2,'-*');
    title('thermal flux');

end

