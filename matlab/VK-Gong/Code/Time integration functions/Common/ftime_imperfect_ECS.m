function [ out ] = ftime_imperfect_ECS( Nphi, Npsi, Ai, H0, H1, H2, C, C1,C2, Tn, e, f_time, rp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ftime_imperfect_ECS Function for the time integration of the Von
%Kármán equation in a imperfect plate either rectangular or circular.

% For rectangular plate, set tnd = 1 and e = Sw^2/rho*E.


%% Energy conserving scheme 




%%
% Pre-conditions:
%      Parameters:
%       -Nphi: Number of transverse modes.
%       - Npsi: Number of Airy functions (in-plane modes)
%      Plate characteristics:
%       - hd: Plate thickness
%       - Ai: Imperfection projection coefficients
%      Non-dimensionalisation variables:
%       - tnd: Non-dimensioning factor
%       - e: Dimensionless stiffness constant e = 12*(1-nu^2)
%      Pre-computed variables: 
%       - H: Matrix of H coefficients. Dimensions (Npsi*Nphi*Nphi).
%           H0: H^i_jk
%           H1: H^i_jk/(\zeta_i)^2
%           H2: H^i_jk/(\zeta_i)^4
%       - mode_t_phys: Frequencies of the transverse modes in radian/second
%               mode_t_phys(:,1) : Mode number
%               mode_t_phys(:,2) : xi
%               mode_t_phys(:,3) : k
%               mode_t_phys(:,4) : n
%               mode_t_phys(:,5) : omega (xi^4=omega^2)
%       - Damping (C, C1, C2): (Nphi*1)
%       - Input force vector: f_time
%       - Output point shape vector: rp
%      Time simulation parameters:
%       - fs: Dimensionless sampling frequency
%       - Tn: Number of iterations
%
%      

tic

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Initialisation   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-linear coefficients%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    H0 = H0(1:Npsi,1:Nphi,1:Nphi)*sqrt(e/2);     % H0 =H
    H1 = H1(1:Npsi,1:Nphi,1:Nphi)*sqrt(e/2);     % H1 =H/zeta^2
    H2 = H2(1:Npsi,1:Nphi,1:Nphi)*sqrt(e/2);     % H0 =H/zeta^4

    
    H0 = reshape(H0,[Npsi*Nphi,Nphi]);
    H1 = reshape(H1,[Npsi*Nphi,Nphi]);
    H2 = reshape(H2,[Npsi*Nphi,Nphi]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Time-Integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

str = sprintf('\n\n-------------------------TIME INTEGRATION-------------------------\n\n');
disp(str);

%initialise q's
q1 = zeros(Nphi,1); %q at timestep n - 1
q2 = zeros(Nphi,1); %q at timestep n - 2
eta1 = zeros(Npsi,1); %eta at timestep n - 1
eta2 = zeros(Npsi,1); %eta at timestep n - 2



%% Output vector initialization (time series)
out = zeros(Tn, size(rp,1));

IterationControl = round(Tn/20);
tLoop = tic;

%main loop
for i=1:Tn
    
    %% Function to be integrated:
    %  qs(n+1)=C\(C1*qs(n)+C2*qs(n-1)-E*Sw^2/(2*rho)*(Gtotal)+ps(n))
    %  The modes have been normalized by norm_modes, we assume that Sw=1.
    
    
    %%% Left side
    
    % (Hmn*Hks)/zeta^4 * (qj(n-1)+aj)*(qk(n-1)+ak)
    t0 = H1*(q1 + Ai);
    t0 = reshape(t0,[Npsi,Nphi]);
    G0 = (t0.'*t0);
    
    mat_imp = C + G0;
    
    %%% Right side
    % (Hmn*Hks)/zeta^4 * qi(n-1)*aj*(qk(n-1)+ak)
    t1 = H1*q1;
    t1 = reshape(t1,[Npsi,Nphi]);
    t4 = t1*Ai;
    Ga = t0.'*t4;
    
    
    % Hks * (qk(n-1)+ak) * (etal(n-2)-etal(n-1))/2
    eta_temp = (eta2-eta1);
    t5 = H0*(q1 + Ai);
    t5 = reshape(t5,[Npsi,Nphi]);
    Gb = (t5.'*eta_temp);
    
    Gtotal = Ga - Gb;
    
    %solve for displacement
    q = mat_imp\(- C1.*q1 - C2.*q2 - Gtotal + f_time(:,i));
    
    %get output
    out(i, :) = rp*q;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Progress control
    if (i == IterationControl)
        Time5 = toc(tLoop);
        EstimatedTime = Time5*20;
        hours = floor(EstimatedTime / 3600);
        EstimatedTime = EstimatedTime - hours * 3600;
        mins = floor(EstimatedTime / 60);
        secs = EstimatedTime - mins * 60;
        disp(['Estimated time for the time integration: ', num2str(hours), ' hours ', num2str(mins), ' minutes ', num2str(secs), ' seconds ']);
    end
    if ( mod(i, IterationControl) == 0) 
        disp(['Progress: ', num2str(i*5/IterationControl), '%']);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Update eta vector
    
    % Hks/zeta^4 * qi(n)*qj(n-1)
    t7 = H2*q;
    t7 = reshape(t7,[Npsi,Nphi]);
    Gc_1 = t7*q1;
    
    % Hks/zeta^4 * aj * (qi(n) + qi(n-1))
    t8 = H2*Ai;
    t8 = reshape(t8,[Npsi,Nphi]);
    Gc_2 = t8*(q + q1);
     
    Gc = (Gc_1 + Gc_2);
    
    eta = -eta1 - Gc;
    
    %update variables
    q2 = q1;
    q1 = q;
    eta2 = eta1;
    eta1 = eta;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
end



str = sprintf('\n\n-------------------END TIME INTEGRATION----------------------------\n\n');
disp(str);
toc







end

