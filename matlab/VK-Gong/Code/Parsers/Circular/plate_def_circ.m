function [Rd, hd, E, BC, e, Nphi, Npsi, scheme, H0, H1, H2, filename, Ai, C, C1, C2, k_t, c_t, xkn, JJ, II, Kkn, rp, tnd, fs, Tsd] = plate_def_circ(PlateCharacteristicsFileName, SimulationParametersFileName, OutputFileName, GammaFileName )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function parses the files indicated in the main script in order to 
% define the plate and create the necessary variables for the code. 

%% Input parameters		
% PlateCharacteristicsFileName: Name of the file containing the plate characteristics parameters. 
% SimulationParametersFileName: Name of the file containing the simulation parameters. 
% OutputFileName: Name used to save the results of the program execution.
% GammaFileName: Name of the file containing the Gamma matrix. If a file with this name does not exist, it will be created during the performance of the code.

%% Output parameters
% Rd: Plate radius in meters.			
% hd: Dimensional plate thickness in meters.
% E: Young modulus in Pa. 
% BC: Type of boundary conditions at the edge. Possible values: 'free', 'clamped', 'elastic'.
% e: Dimensioning or non-dimensioning factor corresponding to $\varepsilon'$.
% Nphi: Number of transverse modes kept in the truncation.
% Npsi: Number of in-plane modes kept in the truncation.
% scheme: Time integration scheme. Possible values: 'ECS', 'verlet'.
% H0: Matrix of size $Npsi\times Nphi\times Nphi$ containing the ${H}$ coupling coefficient tensor, $H(i,p,q)= {H}^i_{pq}$.
% H1: Matrix that contains matrix H0 divided by the eigenfrequency of the corresponding in-plane mode, i.e. $H1(i,p,q)= {H}^i_{pq}/\omega_i = {H}^i_{pq}/\zeta^2_i$.
% H2: Matrix that contains matrix H0 divided by the squared eigenfrequency of the corresponding in-plane mode, i.e. $H2(i,p,q)= {H}^i_{pq}/\omega^2_i= {H}^i_{pq}/\zeta^4_i$.\\	
% filename: String containing the route and the filename to save the results of the code. 
% Ai: Full vector of projection coefficients.
% C: Matrix corresponding to $\left(\frac{1}{k^2}+\frac{C_{ss}}{2k}\right)$ as defined in \cref{eq:4.13}.
% C1: Matrix corresponding to $\left(-\frac{2}{k^2}+K_{ss}\right)$ as defined in \cref{eq:4.14}.
% C2: Matrix corresponding to $\left(\frac{1}{k^2}-\frac{c_s}{k}\right)$.
% k_t: Vector containing the number of nodal diameters of the transverse modes. 
% c_t: Vector containing the configuration of the transverse modes.
% xkn: Vector containing the eigenvalues $\xi$ of the transverse modes.
% JJ: Vector containing $\tilde{J}^f_k(x)$, $\tilde{J}^c_k(x)$ or $\tilde{J}^e_k(x)$ depending on the value of BC. 
% II: Vector containing $\tilde{I}^f_k(x)$, $\tilde{I}^c_k(x)$ or $\tilde{I}^e_k(x)$ depending on the value of BC. 
% Kkn: Vector containing the norm of the transverse modes. 
% rp: Vector containing the modal deformation at the output points. 
% tnd: Non-dimensioning time factor.
% fs: Dimensionless sampling frequency.
% Tsd: Time length of the output signal in seconds.


%% Plate characteristics
if exist(PlateCharacteristicsFileName, 'file') == 2
    load(PlateCharacteristicsFileName);
else
    disp('Plate characteristics file not found. Using preset instead.');
    % Geometrical characteristics
    Rd=0.2; % Plate radius
    hd = .8e-3; % Plate thickness
    
    % Imperfection characteristics
    H = 1e-2; % Imperfection height
    ImperfectionType = 'Spherical'; % Spherical cap
    proj = [];
    modeIndices = [];
    
    tau2 = 0; % Parameter to determine depending on the shape of the cap
    ModeType = 'Axisymmetric'; % Mode type according to its symmetry.
    error_coef = 0.10; % Admitted error from 0 to 1.
    
    % Material parameters
    nu = 0.38;
    E = 2e11;
    rho = 7860;

    % Damping  % X = {'Undamped', 'PowerLaw'}
    X = [];  
    c = zeros(1500, 1);    
    dFac = 0;
    dExp = 0;
    dCons = 0;
    
    % Boundary conditions parameters
    BC = 'elastic'; % Boundary conditions: 'free', 'clamped', 'elastic'
    KR = 10; % Rotational stiffness. KR = Rd*Kr/D;
    KT = 1000; % Translational stiffness. KT = Rd^3*Kt/D;

end

%%%%%%%%%%%%%
% Free equals elastic with KR = 0 and KT = 0
% Clamped equals elastic with KR = inf and KT = inf
switch BC
    case 'free'
        KR = 0;
        KT = 0;
    case 'clamped'
        KR = inf;
        KT = inf;
end

D = E*hd^3/(12*(1-nu^2));
e = 12*(1-nu^2);


%% Simulation parameters
if exist(SimulationParametersFileName, 'file') == 2
    load(SimulationParametersFileName);
else
    disp('Simulation parameters file not found. Using preset instead.');
    
    % Accuracy parameters
    Nphi = 10; %select number of transverse modes
    Npsi = 3; %transverse modes
    NA = 10; %Number of modes considered to compute the eigenfrequencies of the imperfect plate.
    
    % Time simulation parameters
    scheme = 'ECS'; %% Integration scheme: "ECS" or "verlet"
    fsd = 0; %% Sampling frequency
    Tsd = 10; %Simulation time
    
    % Model parameters
    Nr = 400; % Number of discretization points for r
    Nth = 500; % Number of discretization points for \theta
    dr_H = 1e-4; % Integration step used for computing the H coefficients
    
    % Output points
    op = [0.5192 0.8962];
    
    % Damping
    X = [];
    
    c = zeros(Nphi, 1);
    
end

dr=Rd/Nr; % Integration step

Nop =  size(op,1); % Number of output points



%% H coefficients and mode files loading

switch BC
    case 'free'
        modeTfilename = sprintf('Parameters/Mode files/Circular/mode_t_free-nu_%f.mat', nu);
        if exist(modeTfilename, 'file') == 2
            load(modeTfilename);
            mode_t = mode_t_free;
            clear mode_t_free
        else
            dx = 1e-3;
            xmax = 100;
            disp(['Transverse modes file not found. Transverse modes will be computed using dx = ', num2str(dx), 'and xmax = ', num2str(xmax), '. For different values of dx and xmax, execute ComputeTransverseEigenfrequenciesCircular.m externally.' ]);
            [mode_t, ~] = ComputeTransverseEigenfrequenciesCircular( dx, xmax, BC, nu, KR, KT );
        end
        modeLfilename = sprintf('mode_l_free.mat');
        if exist(modeLfilename, 'file') ~= 2
            dx = 1e-3;
            xmax = 100;
            disp(['In-plane modes file not found. In-plane modes will be computed using dx = ', num2str(dx), 'and xmax = ', num2str(xmax), '. For different values of dx and xmax, execute ComputeInPlaneEigenfrequenciesCircular.m externally.' ]);
            ComputeInPlaneEigenfrequenciesCircular( dx, xmax, BC, nu );
        end

        
    case 'clamped'
        modeTfilename = sprintf('mode_t_clamped.mat');
        if exist(modeTfilename, 'file') == 2
            load(modeTfilename);
            mode_t = mode_t_clamped;
            clear mode_t_clamped
        else
            dx = 1e-3;
            xmax = 100;
            disp(['Transverse modes file not found. Transverse modes will be computed using dx = ', num2str(dx), ' and xmax = ', num2str(xmax), '. For different values of dx and xmax, execute ComputeTransverseEigenfrequenciesCircular.m externally.' ]);
            [mode_t, ~] = ComputeTransverseEigenfrequenciesCircular( dx, xmax, BC, nu, KR, KT );
        end
        
        modeLfilename = sprintf('mode_l_clamped-nu_%d.mat', nu);
        if exist(modeLfilename, 'file') ~= 2
            dx = 1e-3;
            xmax = 100;
            disp(['In-plane modes file not found. In-plane modes will be computed using dx = ', num2str(dx), ' and xmax = ', num2str(xmax), '. For different values of dx and xmax, execute ComputeInPlaneEigenfrequenciesCircular.m externally.' ]);
            ComputeInPlaneEigenfrequenciesCircular( dx, xmax, BC, nu );
        end
        
    case 'elastic'
        modeTfilename = sprintf('Parameters/Mode files/Circular/mode_t_elastic-nu_%f-KR_%f-KT_%f.mat', nu, KR, KT);
        if exist(modeTfilename, 'file') == 2
            load(modeTfilename);
            mode_t = mode_t_elastic;
            clear mode_t_elastic
        else
            dx = 1e-3;
            xmax = 100;
            disp(['Transverse modes file not found. Transverse modes will be computed using dx = ', num2str(dx), ' and xmax = ', num2str(xmax), '. For different values of dx and xmax, execute ComputeTransverseEigenfrequenciesCircular.m externally.' ]);
            [mode_t, ~] = ComputeTransverseEigenfrequenciesCircular( dx, xmax, BC, nu, KR, KT );
        end
        
        modeLfilename = sprintf('mode_l_elastic.mat');
        if exist(modeLfilename, 'file') ~= 2
            dx = 1e-3;
            xmax = 100;
            disp(['In-plane modes file not found. In-plane modes will be computed using dx = ', num2str(dx), ' and xmax = ', num2str(xmax), '. For different values of dx and xmax, execute ComputeInPlaneEigenfrequenciesCircular.m externally.' ]);
            ComputeInPlaneEigenfrequenciesCircular( dx, xmax, BC, nu );
        end
        
    otherwise
        disp('Unknown boundary conditions');
end

disp('Mode files loaded');

[H0, H1, H2] = LoadHTensorCircular(BC, Nphi, Npsi, nu, KR, KT, dr_H, scheme);

disp('H tensors loaded');

%% Output file parameters
%OutputFileName = 'Cymb_Test1';
mkdir(OutputFileName);
filename = sprintf('%s/%s_%s_%s', OutputFileName,OutputFileName, scheme, BC);

%% Imperfection modelling
if isempty(proj)
    if H == 0
        Ai = zeros(Nphi,1); % Projection coefficients
        proj = [];
        modeIndices = [];
        R = inf;
    else

        [ ~, ~, ~, proj, modeIndices, ~, R] = ProjectionOfTheImperfectionCircular( H,hd,Rd, ImperfectionType, tau2, Nr,Nth,nu,KR,error_coef,ModeType, mode_t );
        Ai = zeros(Nphi,1); % Projection coefficients
        Ai(modeIndices)=proj;
        Ai=Ai(1:Nphi); %%In case that Nphi<max(modeIndices)
        
        disp('Performed computation of projection coefficients');
        
    end
save(sprintf('%s-Imperfection_Parameters.mat',filename),'Ai','modeIndices', 'proj','R');    
end


%% Dimensionless eigenfrequencies

om = mode_t(:,2).^2;    % Eigenfrequencies


%% Computation of eigenfrequencies of the imperfect plate
if NA ~= 0
    GammaFileName = sprintf('Parameters/H files/Circular/%s',GammaFileName);
    if exist(GammaFileName, 'file') ~= 2
        disp('Gamma file not found. Creating it');
        GammaTensor( H1, GammaFileName, NA, Npsi );
        
    end
    
    
    [Omega , ~] = ComputeEigenfrequenciesImperfectPlate(proj,NA, modeIndices, om, e, GammaFileName );

    disp('Computed eigenfrequencies of the imperfect plate');
    
    Chi=(12*(1-nu^2)*Rd^4)/(R^2*hd^2); % Equivalent curvature, only useful for spherical imperfections.

    save (sprintf('%s-Eigenfrequencies-Curvature.mat',filename),'Omega','Chi');
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Nondimensioning %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


R = Rd/Rd;        % Dimenssionless radius

dr = dr/Rd;       % Dimensionless integration step

tnd = Rd^2*sqrt(rho*hd/D);   % Nondimensionalization time coefficient t0

om_dim = om/tnd; % Dimensional omega
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% sampling rate (Hz),
fsd_lim = om_dim(Nphi); %limit frequency

if (fsd < fsd_lim)
    disp(['Warning: Input sampling rate is below the stability limit. New fsd = ', num2str(fsd_lim)]);
    fsd = fsd_lim;
end

fs = fsd*tnd;     % Dimensionless sampling rate
fs = round(fs); %sampling rate


%% damping ratios
if ~isempty(X)
    [c] = c_preset(X,om_dim(1:Nphi),Nphi, dFac, dExp, dCons);
end

% Dimensionless damping
c=Rd^2/sqrt(rho*hd*D)*c*rho*hd;


%timestep
k = 1/fs;






%% Initialise constant matrices for main loop
C = (1/k^2 + c(1:Nphi)/k/(1));
C1 = zeros(Nphi,1);
C2 = zeros(Nphi,1);

switch scheme 
    case 'ECS' %% Energy Conserving Scheme
        C = diag(C);

        C1 = (-2/k^2 + (om(1:Nphi)/(1)).^2);

        C2 = (1/k^2 - c(1:Nphi)/k/(1));

    case 'verlet' % Stormer - verlet
        C1 = (-2/k^2 + (om(1:Nphi)/(1)).^2)./C;

        C2 = (1/k^2 - c(1:Nphi)/k/(1))./C;

    otherwise
        disp('Unknown scheme');
end

switch BC
    case {'free', 'elastic'}
        %% Compute rp factors applied to input 
        rp = zeros(Nop,Nphi);
        J0 = zeros(1,Nphi);
        J1 = zeros(1,Nphi);
        J2 = zeros(1,Nphi);
        I0 = zeros(1,Nphi);
        I1 = zeros(1,Nphi);
        I2 = zeros(1,Nphi);
        JJ0 = zeros(Nop,Nphi);
        II0 = zeros(Nop,Nphi);
        Kkn = zeros(1,Nphi);
        Jtild = zeros(1,Nphi);
        Itild = zeros(1,Nphi);


        c_t= mode_t(1:Nphi,5);
        k_t = mode_t(1:Nphi,3)'; 
        xkn = sqrt(om)';
        clear mode_t


        for ii = 1 : Nphi

            J0(ii)=besselj(k_t(ii),xkn(ii));
            J1(ii)=besselj(k_t(ii)-1,xkn(ii));
            J2(ii)=besselj(k_t(ii)-2,xkn(ii));

            I0(ii)=besseli(k_t(ii),xkn(ii)); 
            I1(ii)=besseli(k_t(ii)-1,xkn(ii));
            I2(ii)=besseli(k_t(ii)-2,xkn(ii));

            Jtild(ii)=xkn(ii)^2*J2(ii)+((nu-2*k_t(ii)+1)*xkn(ii) + KR)*J1(ii)+(k_t(ii)*(k_t(ii)+1)*(1-nu)-KR*k_t(ii))*J0(ii);
            Itild(ii)=xkn(ii)^2*I2(ii)+((nu-2*k_t(ii)+1)*xkn(ii) + KR)*I1(ii)+(k_t(ii)*(k_t(ii)+1)*(1-nu)-KR*k_t(ii))*I0(ii);

            JJ0(:,ii) = besselj(k_t(ii),xkn(ii)*op(:,2));
            II0(:,ii) = besseli(k_t(ii),xkn(ii)*op(:,2));

            rp(:,ii) = (JJ0(:,ii) - (Jtild(ii)*II0(:,ii)/(Itild(ii)))).*cos(k_t(ii)*op(:,1)-(c_t(ii)-1)/2*pi);  % Modal deformation at output point  

            % Normalisation 
            Kkn(ii) = norm_modes(k_t(ii),xkn(ii),R,dr,nu, KR, BC);

            rp(:,ii) = rp(:,ii) * Kkn(ii);
            
        end
        
        JJ = Jtild;
        II = Itild;
        
    case 'clamped'
        rp = zeros(Nop,Nphi);
        Jkn = zeros(1,Nphi);
        Ikn = zeros(1,Nphi);
        Kkn = zeros(1,Nphi);       
        
        c_t= mode_t(1:Nphi,5);
        k_t = mode_t(1:Nphi,3)';
        xkn = sqrt(om)';
        clear mode_t      

        for ii = 1 : Nphi

            Jkn(ii)=besselj(k_t(ii),xkn(ii));
            Ikn(ii)=besseli(k_t(ii),xkn(ii));

            J=besselj(k_t(ii),xkn(ii)*op(:,2));
            I=besseli(k_t(ii),xkn(ii)*op(:,2));

            rp(:,ii) = (Ikn(ii)*J-Jkn(ii)*I).*cos(k_t(ii)*op(:,1)-(c_t(ii)-1)/2*pi); % Modal deformation at output point
           

            % Normalisation 
            Kkn(ii) = norm_modes(k_t(ii),xkn(ii),R,dr,nu, KR, BC);

            rp(:,ii) = rp(:,ii) * Kkn(ii);


        end
        
        JJ = Jkn;
        II = Ikn;
    
    otherwise
        disp('Unknown boundary conditions');

end
        


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp ('Plate input files parsed');

end

