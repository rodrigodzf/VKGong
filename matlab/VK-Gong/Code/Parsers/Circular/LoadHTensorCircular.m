function [H0, H1, H2] = LoadHTensorCircular(BC, Nphi, Npsi, nu, KR, KT, dr_H, scheme)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that loads the H file or creates it if necessary. 

list = dir('Parameters/H files/Circular');
Found = 0;

for i = 1:length(list)
    Name = list(i).name;
    switch BC
        case 'free'
            if regexp(Name, 'H_free')
                [Nphif, Npsif, nuf, dr_Hf] = strread(Name, 'H_free-Nphi_%d-Npsi_%d-nu_%f-dr_%f.mat'); 
                if ((Nphif>=Nphi)&&(Npsif>=Npsi)&&(nuf==nu)&&(dr_Hf<=dr_H))
                    switch (scheme)
                        case 'ECS'
                            load(Name);

                        case 'verlet'
                            load(Name, 'H1');
                            H0 = [];
                            H2 = [];
 
                    end
                    Found = 1;
                    break;
                end
            end

        case 'clamped'
            if regexp(Name, 'H_clamped')
                [Nphif, Npsif, nuf, dr_Hf] = strread(Name, 'H_clamped-Nphi_%d-Npsi_%d-nu_%f-dr_%f.mat');
                if ((Nphif>=Nphi)&&(Npsif>=Npsi)&&(nuf==nu)&&(dr_Hf<=dr_H))
                    switch (scheme)
                        case 'ECS'
                            load(Name);

                        case 'verlet'
                            load(Name, 'H1');
                            H0 = [];
                            H2 = [];

                    end
                    Found = 1;
                    break;
                end
            end


        case 'elastic'
            if regexp(Name, 'H_elastic')
                [Nphif, Npsif, nuf, dr_Hf, KRf, KTf] = strread(Name, 'H_elastic-Nphi_%d-Npsi_%d-nu_%f-dr_%f-KR_%f-KT_%f.mat');
                if ((Nphif>=Nphi)&&(Npsif>=Npsi)&&(nuf==nu)&&(dr_Hf<=dr_H)&&(KRf==KR)&&(KTf==KT))
                    switch (scheme)
                        case 'ECS'
                            load(Name);

                        case 'verlet'
                            load(Name, 'H1');
                            H0 = [];
                            H2 = [];

                    end
                    Found = 1;
                    break;
                end
            end

        otherwise
            disp('Unknown boundary conditions');
    end
    if Found == 1
        break;
    end
end

if ~Found
    disp('H file not found. Computing H coefficients.');
    disp('Warning: This may take a long time.');   
    
    switch (scheme)
        case 'ECS'
            [ H0, H1, H2 ] = H_tensorCircular( Nphi, Npsi, BC, nu, KR, KT, dr_H );
            H0 = H0(1:Npsi, 1:Nphi, 1:Nphi);
            H1 = H1(1:Npsi, 1:Nphi, 1:Nphi);
            H2 = H2(1:Npsi, 1:Nphi, 1:Nphi); 

        case 'verlet'
            [ ~, H1, ~ ] = H_tensorCircular( Nphi, Npsi, BC, nu, KR, KT, dr_H );
            H1 = H1(1:Npsi, 1:Nphi, 1:Nphi);
            H0 = [];
            H2 = [];

    end
    
    
else
    switch (scheme)
        case 'ECS'
            H0 = H0(1:Npsi, 1:Nphi, 1:Nphi);
            H1 = H1(1:Npsi, 1:Nphi, 1:Nphi);
            H2 = H2(1:Npsi, 1:Nphi, 1:Nphi); 
            
        case 'verlet'
            H1 = H1(1:Npsi, 1:Nphi, 1:Nphi);
            H0 = [];
            H2 = [];
            
    end
end



end

