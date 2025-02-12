function [ U,V,Imperfection,proj,modeIndices,Approximation,Rc ] = ProjectionOfTheImperfectionCircular( H,hd,Rd, ImperfectionType, tau2, Nr,Nth,nu,KR,error_coef,ModeType, mode_t )
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

%% Input parameters
% H: Cap height
% hd: Thickness
% Rd: Cap radius
% ImperfectionType:
%   Spherical
%   Parabolic 
% tau2: Parabola order
% nbpointsr: Number of points of discretization in r
% nbpointsth: Number of points of discretization in th
% nu: Poisson ratio.
% KR: Normalized rotational stiffness with respect to bending stiffness. KR = Kr/D  Only used when BC = 'elastic'.
% error_coef: Top error admitted in the approximation of the imperfection
% ModeType: Type of modes considered in the approximation of the imperfection. Possible values: 'All', 'Axisymmetric'.
% mode_t: Vector containing the information corresponding to the translational vibration modes.
%      <1:Index> <2:x_i> <3:k> <4:n> <5:c> <6:x^2_i>

%% Output parameters:
% U,V: Coordinates of points in Imperfection
% Imperfection: Profile surface
% proj: Projection coefficients
% modeIndices: Indices of modes used to approximate the imperfection.
% Approximation: Profile obtained using the projection coefficients. 
% Rc: Curvature radius of the imperfection. 


Rc = (H^2 + Rd^2)/(2*H);

[U,V,Imperfection]=AxisymmetricCap(H,Rd,ImperfectionType,Nr,Nth,tau2); % Computation of the imperfection profile in polar coordinates. 

[proj, modeIndices,Approximation, ~]=ComputationOfTheProjectionCoefficientsCircular(Imperfection/hd,error_coef,nu,KR,ModeType, mode_t); 



end
