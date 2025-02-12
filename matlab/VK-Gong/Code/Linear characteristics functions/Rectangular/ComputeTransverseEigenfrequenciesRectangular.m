function [mode_t] = ComputeTransverseEigenfrequenciesRectangular( BC, Lx, Ly,  Nphi  )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the transverse eigenvalues \xi (\xi^2 = \omega) 
% for a corresponding set of boundary conditions in a rectangular plate. 

N = ceil(sqrt(Nphi))+1;

mode_t = zeros(N^2, 3);


switch BC            
    case {'SimplySupported'}
        filename = sprintf('Parameters/Mode files/Rectangular/mode_t_SS-Lx_%f-Ly_%f-Nphi_%d.mat', Lx, Ly, Nphi);
        if exist(filename, 'file') == 2
            load(filename);
        else
            [kx, ky] = meshgrid(1:Nphi, 1:Nphi);

            kx = reshape(kx, Nphi^2, 1);
            ky = reshape(ky, Nphi^2, 1);

            gf = ((kx*pi/Lx).^2 + (ky*pi/Ly).^2);

            mode_t = [kx ky gf];

            mode_t = sortrows(mode_t,3);

            mode_t = [(1:size(mode_t,1))' mode_t];

            mode_t = mode_t(1:Nphi, :);

            save(filename, 'mode_t');
        end

    otherwise
            disp('Unknown boundary conditions');

end
        

    
end

