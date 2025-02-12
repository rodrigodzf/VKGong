function [ out ] = ftime_imperfect_verlet( Nphi, Npsi, Ai, H1, C, C1, C2, Tn, e, f_time, rp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ftime_imperfect_stable Function for the time integration of the Von
%Kármán equation in a imperfect plate either rectangular or circular.

% For rectangular plate, set tnd = 1 and e = Sw^2/rho*E.


%% Stormer-Verlet scheme (Unstable)

% Dimless stands for dimensionless.
% Modular is used to distinguish from the version that does not separate
% plate_def and score. 


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
%           H1: H^i_jk/(\zeta_i)^2
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-linear coefficients%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    H1 = H1(1:Npsi,1:Nphi,1:Nphi)*sqrt(e/2);     % H1 =H/zeta^2

    H1 = reshape(H1,[Npsi*Nphi,Nphi]);
    


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Time-Integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

str = sprintf('\n\n-------------------------TIME INTEGRATION-------------------------\n\n');
disp(str);

%initialise q's
q1 = zeros(Nphi,1); %q at timestep n - 1
q2 = zeros(Nphi,1); %q at timestep n - 2



%% Output vector initialization (time series)
out = zeros(Tn, size(rp,1));

IterationControl = round(Tn/20);
tLoop = tic;

%main loop
for i=1:Tn
    
   
    % (Hmn*Hks)/2zeta^4 * (qk(n)+ak)*qm(n)*(qn(n)+2*an) = e/2 *(qk(n)+ak(n))*qm(n)*(qn(n)+2*an)
    t0 = H1*(q1 + Ai);
    t0 = reshape(t0,[Npsi,Nphi]);
    t1 = H1*q1;
    t1 = reshape(t1,[Npsi,Nphi]);
    t2 = t1*(q1 + 2*Ai);
    G = t0.'*t2;
    
    %solve for displacement
    q = - C1.*q1 - C2.*q2 - G./C + f_time(:,i);
    
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
    
    %update variables
    q2 = q1;
    q1 = q;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
end



str = sprintf('\n\n-------------------END TIME INTEGRATION----------------------------\n\n');
disp(str);
toc







end

