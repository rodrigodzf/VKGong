function [U,V,Imperfection]=AxisymmetricCap(H,R,ImperfectionType,Nr,Nth, tau2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Building of the imperfection profile as a 3D surface.
%% Input parameters
% H: Height of the imperfection.
% R: Plate radius. 
% Nr: Number of discretization points in r
% Nth: Number of discretization points in theta discr�tisation en theta
% tau2: Curve ratio for the imperfections with two different curves or the types "chinessecymbal" / Parabola order
% ImperfectionType:
%   Spherical
%   Parabolic 



%% The profile of the imperfection is computed in cartesian coordinates and after transformed to polar coordinates. 

discrtheta=Nth;
r=0:R/(Nr-1):R;
theta=0:2*pi/(discrtheta-1):2*pi;

switch ImperfectionType
    case 'Spherical'
        % According to the definitions included in the thesis of Cedric
        % Camier and the article doi:10.1016/j.euromechsol.2008.11.005.  
        %It should be normalized by the thickness. 
        Rc = (H^2 + R^2)/(2*H);
        y=-Rc+sqrt(Rc^2-r.^2);
        y0 = pi*((2/3)*(Rc^3-(Rc^2-R^2)^(3/2))-Rc*R^2); %% Closed form of the integral. 
        y=y-y0/(pi*R^2);

        
    case 'Parabolic'
        y=-H*r.^tau2;
        y0 = -2*pi*H*R^(tau2 + 2)/(tau2 + 2);
        y=y-y0/(pi*R^2);
    

end
    
Imperfection=y'*ones(size(theta));
U=r'*cos(theta);
V=r'*sin(theta);
% figure
% surf(U,V, Imperfection)
% shading interp
end

