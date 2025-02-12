%% This main file is used to execute the simulation functions for the nonlinear vibrations of thin plates. 

%% Input parameters files
PlateCharacteristicsFileName = 'PlateCharacteristicsCircular.mat'; % Physical characteristics of the plate: Dimensions, imperfection profile, material and boundary conditions.  
SimulationParametersFileName = 'SimulationParametersCircular.mat'; % Parameters related to the simulation: Time length, scheme, number of modes, output points, accuracy.
GammaFileName = 'GammaCircular-Nphi_NPHICircular.mat'; % Name of the file containing the Gamma Tensor.  
ScoreFileName = 'ScoreParametersCircular.mat'; % Characteristics of the excitation. 
OutputFileName = 'ResultsCircular'; % Name of the results files and folder. 

%% Simulation setup 
[Rd, hd, E, BC, e, Nphi, Npsi, scheme, H0, H1, H2, filename, Ai, C, C1, C2, k_t, c_t, xkn, JJ, II, Kkn, rp, tad, fs, Tsd] = plate_def_circ(PlateCharacteristicsFileName, SimulationParametersFileName, OutputFileName, GammaFileName );

[ f_time, Tn ] = score_circ(ScoreFileName, Rd, hd, E, BC, e, Nphi, scheme, C, k_t, c_t, xkn, JJ, II, Kkn, tad, fs, Tsd);

%% Time simulation
switch scheme
    case 'ECS'
        [ out_nd ] = ftime_imperfect_ECS( Nphi, Npsi, Ai, H0, H1, H2, C, C1,C2, Tn, e, f_time, rp);
        
    case 'verlet'
        [ out_nd ] = ftime_imperfect_verlet( Nphi, Npsi, Ai, H1, C, C1, C2, Tn, e, f_time, rp);
    
    otherwise
        disp('Unknown scheme');
end

% Save results
fsd = round(fs/tad); % Dimensioned sampling frequency
out = out_nd*hd; % Dimensioned output displacement
out_vel = diff(out,1,1)*fsd; % Dimensioned output velocity
save(sprintf('%s.mat',filename), 'out', 'out_vel',  'fsd');

%% Generate audio from output (Comment this block if the code is not used for sound synthesis purposes)
for i = 1:size(rp,1)
    audiowrite(sprintf('%s-Op_%d.wav',filename, i),out_vel(:,i)/1.1/max(abs(out_vel(:,i))),fsd);
end
