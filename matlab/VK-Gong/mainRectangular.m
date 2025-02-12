%% This main file is used to execute the simulation functions for the nonlinear vibrations of thin plates. 

%% Input parameters files
PlateCharacteristicsFileName = 'PlateCharacteristicsRectangular.mat'; % Physical characteristics of the plate: Dimensions, imperfection profile, material and boundary conditions.  
SimulationParametersFileName = 'SimulationParametersRectangular.mat'; % Parameters related to the simulation: Time length, scheme, number of modes, output points, accuracy.
GammaFileName = 'GammaRectangular.mat'; % Name of the file containing the Gamma Tensor.  
ScoreFileName = 'ScoreParametersRectangular.mat'; % Characteristics of the excitation. 
OutputFileName = 'ResultsRectangular'; % Name of the results files and folder. 

%% Simulation setup 
[Lx, Ly, hd, E, rho, BC, e,  Nphi, Npsi, scheme, H0, H1, H2, filename, Ai, C, C1, C2, kx, ky, om_dim, rp, tad, fsd, Tsd] =plate_def_rect(PlateCharacteristicsFileName, SimulationParametersFileName, OutputFileName, GammaFileName );

[ f_time, Tn ] = score_rect( ScoreFileName, Lx, Ly, hd, rho, kx, ky, BC, Nphi, scheme, C, fsd, Tsd);

%% Time simulation
switch scheme
    case 'ECS'
        [ out ] = ftime_imperfect_ECS( Nphi, Npsi, Ai, H0, H1, H2, C, C1,C2, Tn, e, f_time, rp);
        
    case 'verlet'
        [ out ] = ftime_imperfect_verlet( Nphi, Npsi, Ai, H1, C, C1, C2, Tn, e, f_time, rp);
    
    otherwise
        disp('Unknown scheme');
end

%% Save results
out_vel = diff(out,1,1)*fsd; % Output velocity
save(sprintf('%s.mat',filename), 'out', 'out_vel',  'fsd');

%% Generate audio from output (Comment this block if the code is not used for sound synthesis purposes)
for i = 1:size(rp,1)
    audiowrite(sprintf('%s-Op_%d.wav',filename, i),out_vel(:,i)/1.1/max(abs(out_vel(:,i))),fsd);
end
