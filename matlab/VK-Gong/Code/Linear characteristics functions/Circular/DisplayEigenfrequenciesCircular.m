function [ om_dim, f_dim, om_ndim ] = DisplayEigenfrequenciesCircular( PlateCharacteristicsFileName, SimulationParametersFileName, GammaFileName )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters loading
if exist(PlateCharacteristicsFileName, 'file') == 2
    load(PlateCharacteristicsFileName);
else
    disp('Plate Characteristics file not found');
    return
end

if exist(SimulationParametersFileName, 'file') == 2
    load(SimulationParametersFileName);
else
    disp('Simulation Parameters file not found');
    return
end

switch BC
    case 'free'
        modeTfilename = sprintf('mode_t_free-nu_%f.mat', nu);
        load(modeTfilename);
        mode_t = mode_t_free;
        clear mode_t_free
        
    case 'clamped'
        load mode_t_clamped
        mode_t = mode_t_clamped;
        clear mode_t_clamped  
        
    case 'elastic'
        modeTfilename = sprintf('mode_t_elastic-nu_%f-KR_%f-KT_%f.mat', nu, KR, KT);
        load(modeTfilename);
        mode_t = mode_t_elastic;
        clear mode_t_elastic

    otherwise
        disp('Unknown boundary conditions');
end

%% Intermediate parameters
om = mode_t(:,2).^2;    % Eigenfrequencies
D = E*hd^3/(12*(1-nu^2));
e = 12*(1-nu^2);
tnd = Rd^2*sqrt(rho*hd/D);   % Nondimensionalization time coefficient t0


%% Imperfection modelling
if isempty(proj)
    if H == 0
        proj = [];
        modeIndices = [];
        
    else

        [ ~, ~, ~, proj, modeIndices, ~, ~] = ProjectionOfTheImperfectionCircular( H,hd,Rd, ImperfectionType, tau2, Nr,Nth,nu,KR,error_coef,ModeType, mode_t );

    end
     
end


%% Computation of eigenfrequencies

if H ~= 0
    GammaFileName = sprintf('Parameters/H files/Circular/%s',GammaFileName);
    if exist(GammaFileName, 'file') ~= 2
        disp('Gamma file not found. Creating it');
        
        [~, H1, ~] = LoadHTensorCircular(BC, Nphi, Npsi, nu, KR, KT, dr_H, scheme);
        GammaTensor( H1, GammaFileName, NA, Npsi );
        
    end
        
    
    [om_ndim , ~] = ComputeEigenfrequenciesImperfectPlate(proj,NA, modeIndices, om, e, GammaFileName );
else
    
    om_ndim = om;
end

    om_dim = om_ndim/tnd; % Dimensional omega
    
    f_dim = om_dim / 2 / pi; % Dimensional frequency

end

