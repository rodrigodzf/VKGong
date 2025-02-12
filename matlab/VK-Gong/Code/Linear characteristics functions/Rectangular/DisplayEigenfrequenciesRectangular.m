function [ om_dim, f_dim ] = DisplayEigenfrequenciesRectangular( PlateCharacteristicsFileName, SimulationParametersFileName, GammaFileName )
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

[mode_t] = ComputeTransverseEigenfrequenciesRectangular( BC, Lx, Ly,  Nphi  );

%% Intermediate parameters

D = E*hd^3/(12*(1-nu^2));
e = (Lx*Ly/4)/rho*E;
om = sqrt(D/rho/hd)*mode_t(:, 4);


%% Imperfection modelling
if isempty(proj)
    if H == 0       
        proj = [];
        modeIndices = [];

    else

        [ ~,~,~,proj,modeIndices,~ ]  = ProjectionOfTheImperfectionRectangular( H, Lx, Ly, Nx, Ny, BC, error_coef, ModeType, ImperfectionType, xWidth, yWidth, Nphi);
        
    end
   
end
%% Computation of eigenfrequencies

if H ~= 0
    GammaFileName = sprintf('Parameters/H files/Circular/%s',GammaFileName);
    if exist(GammaFileName, 'file') ~= 2
        disp('Gamma file not found. Creating it');
        
        [~, H1, ~] =  LoadHTensorRectangular(BC, Nphi, Npsi, Lx, Ly, mode_t);
        GammaTensor( H1, GammaFileName, NA, Npsi );
        
    end
        
    
    [om_dim , ~] = ComputeEigenfrequenciesImperfectPlate(proj,NA, modeIndices, om, e, GammaFileName );
else
    
    om_dim = om;
end
    
    
    f_dim = om_dim / 2 / pi; % Dimensional frequency

end

