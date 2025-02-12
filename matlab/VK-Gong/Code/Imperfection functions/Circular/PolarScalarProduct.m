function [ I ] = PolarScalarProduct( f1,f2,r_min,r_max,theta_min,theta_max )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function that calculates the scalar product between two functions f1 and
% f2 in polar coordinates. The surface domain of integration is set by the
% boundaries r_min, r_max,theta_min and theta_max. 


r=linspace(r_min,r_max,size(f1,1))';
theta=linspace(theta_min,theta_max,size(f1,2));


r_conversion=repmat(r,1,length(theta));



I=trapz(theta,trapz(r,f1.*f2.*r_conversion,1));



end

