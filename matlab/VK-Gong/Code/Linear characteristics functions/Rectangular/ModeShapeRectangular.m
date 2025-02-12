function [ X, Y, phi ] = ModeShapeRectangular( BC, kx, ky, Lx, Ly, Nx, Ny )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the transverse mode shapes of rectangular plates with boundary
% conditions according to BC. 

switch BC
    case 'SimplySupported'
        x = linspace(0, Lx, Nx);
        y = linspace(0, Ly, Ny);

        [X, Y] = meshgrid(x, y);

        phi = sin(kx*pi*X/Lx).*sin(ky*pi*Y/Ly);

        phi = phi';
        
    otherwise
        disp('Unknown boundary conditions');
end

end

