function [ U,V,Imperfection,proj,modeIndices,Approximation ] = ProjectionOfTheImperfectionRectangular( H, Lx, Ly, Nx, Ny, BC, error_coef, ModeType, ImperfectionType, xWidth, yWidth, Nphi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ProjectionOfTheImperfection This function computes the coefficients to
%project the imperfection of the plate on the modes of the transverse
%motion. 

%% Input parameters:
% H: Deformation height
% Lx, Ly: Plate dimensions
% Nx: Number of points of discretization in x
% Ny: Number of points of discretization in y
% error_coef: Top error admitted in the approximation of the imperfection
% ModeType: Types of modes defined as follows
%   'All'
% ImperfectionType: Type of imperfection
%   '2DRaisedCosine'
% Nphi: Number of transverse modes.

%% Output parameters:
% U,V: Coordinates of points in Imperfection
% Imperfection: Profile surface
% proj: Projection coefficients
% modeIndices: Indices of modes used to approximate the imperfection.
% Approximation: Profile obtained using the projection coefficients. 


[U, V, Imperfection] = RectangularImperfection( Lx, Ly, H, Nx, Ny, ImperfectionType, xWidth, yWidth); % Computation of the imperfection profile

[proj, modeIndices, Approximation, ~] = ComputationOfTheProjectionCoefficientsRectangular( Imperfection, error_coef, ModeType, BC, Nphi, Lx, Ly );



end
