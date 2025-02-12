function [U,V,Imperfection]=RectangularImperfection(Lx,Ly,H, Nx, Ny, ImperfectionType, xWidth, yWidth)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Building of the imperfection profile as a 3D surface. 

%% Input parameters:
% H: Deformation height
% Lx, Ly: Plate dimensions
% Nx: Number of points of discretization in x
% Ny: Number of points of discretization in y
% ImperfectionType: Type of imperfection
%   '2DRaisedCosine'


%% Output parameters:
% U,V: Coordinates of points in Imperfection
% Imperfection: Profile surface


x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);

Imperfection = zeros(Nx, Ny);

[U, V] = meshgrid(x,y);

U = U';
V = V';


switch ImperfectionType
    case '2DRaisedCosine'
        % x Direction
        xVec = zeros(Nx,1);
        X0 = round(Nx/2);
        fsx = Nx/Lx;
        n =  (floor(X0 - fsx*xWidth) : floor(X0 + fsx*xWidth))';
        xVec_temp = 1/2*(1 + cos(pi*((n-1)/fsx/xWidth - X0/fsx/xWidth)));
        xVec(floor(X0 - fsx*xWidth +1) : floor(X0 + fsx*xWidth+1)) = xVec_temp;
        
        % y Direction
        yVec = zeros(Ny,1);
        Y0 = round(Ny/2);
        fsy = Ny/Ly;
        n =  (floor(Y0 - fsy*yWidth) : floor(Y0 + fsy*yWidth))';
        yVec_temp = 1/2*(1 + cos(pi*((n-1)/fsy/yWidth - Y0/fsy/yWidth)));
        yVec(floor(Y0 - fsy*yWidth+1) : floor(Y0 + fsy*yWidth+1)) = yVec_temp;

        Imperfection = xVec(1:Nx)*yVec(1:Ny)';
        %Imperfection = Imperfection - CartesianScalarProduct(Imperfection,ones(size(Imperfection)),0,Lx,0,Ly)/(Lx*Ly);         
        Imperfection = H*Imperfection;
        
end

% figure
% surf(U, V, Imperfection);
% shading interp;


end

