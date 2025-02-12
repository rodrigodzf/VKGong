function [H0, H1, H2] = LoadHTensorRectangular(BC, Nphi, Npsi, Lx, Ly, mode_t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loads the H file and if nonexistent, computes the H matrices and saves
% them. 

list = dir('Parameters/H files/Rectangular');
Found = 0;

for i = 1:length(list)
    Name = list(i).name;
    switch BC
        case 'SimplySupported'
            if regexp(Name, 'H_SS')
                [Nphif, Npsif, Lxf, Lyf] = strread(Name, 'H_SS-Nphi_%d-Npsi_%d-Lx_%f-Ly_%f.mat'); 
                if ((Nphif>=Nphi)&&(Npsif>=Npsi)&&(Lxf==Lx)&&(Lyf==Ly))
                    load(Name);
                    Found = 1;
                    break;
                end
            end

        otherwise
            disp('Unknown boundary conditions');
    end
    if Found == 1
        H0 = reshape(H0,[Npsi*Nphi,Nphi]);
        H1 = reshape(H1,[Npsi*Nphi,Nphi]);
        H2 = reshape(H2,[Npsi*Nphi,Nphi]);
        break;
    end
end

if ~Found
    disp('H file not found. Computing H coefficients.');
    [coeff0, coeff1, coeff2] = AiryStressFactorsCalculation(BC, Npsi,Lx,Ly);
    [ H0, H1, H2 ] = H_tensorRectangular( coeff0, coeff1, coeff2, Nphi, Npsi, Lx, Ly, mode_t, BC);
    
    
    H0 = H0(1:Npsi, 1:Nphi*Nphi);
    H1 = H1(1:Npsi, 1:Nphi*Nphi);
    H2 = H2(1:Npsi, 1:Nphi*Nphi);   

    H0 = reshape(H0,[Npsi, Nphi, Nphi]);
    H1 = reshape(H1,[Npsi, Nphi, Nphi]);
    H2 = reshape(H2,[Npsi, Nphi, Nphi]);
    
    switch BC
        case 'SimplySupported'
            save(sprintf('Parameters/H files/Rectangular/H_SS-Nphi_%d-Npsi_%d-Lx_%f-Ly_%f.mat', Nphi, size(H0,1), Lx, Ly), 'H0', 'H1', 'H2');
            
        otherwise
        disp('Unknown boundary conditions');
    end       
end


end

