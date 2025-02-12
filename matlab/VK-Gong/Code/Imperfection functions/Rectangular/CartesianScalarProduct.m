function [ I ] = CartesianScalarProduct( f1,f2,x_min,x_max,y_min,y_max )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that calculates the scalar product between two functions f1 and
% f2 in cartesian coordinates. The surface domain of integration is set by the
% boundaries x_min, x_max,y_min and y_max. 

x=linspace(x_min,x_max,size(f1,1))';
y=linspace(y_min,y_max,size(f1,2));

I=trapz(y,trapz(x,f1.*f2,1)); 



end

