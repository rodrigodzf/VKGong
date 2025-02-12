function [ zer ] = FindZeros( x, f )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function is used to find the position of the zeros of function f in
% corresponding to vector x

%% Input parameters
% x: Vector of abscissa points
% f: Function to be evaluated

%% Output parameters
% zer: Vector containing the values of x where f is null

fdec=f(2:end);
fdec(end+1)=0;

indice=find((sign(f)+sign(fdec)==0) & (sign(f)~=0)); % Finds the points where the function changes its sign


indiceinc=indice+ones(size(indice));

zer=(x(indice).*f(indiceinc)-x(indiceinc).*f(indice))./(f(indiceinc)-f(indice));


end

