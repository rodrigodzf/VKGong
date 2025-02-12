function [proj, modeIndices, Approximation, error] = ComputationOfTheProjectionCoefficientsRectangular( Imperfection, error_coef, ModeType, BC, Nphi, Lx, Ly )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the coefficients of the series to expansion to 
% express the imperfection profile in terms of the modal basis.


%% Input parameters:
    % Imperfection: Imperfection profile discretized in cartesian coordinates.
    %               size(Imperfection)=[nbpointsx, nbpointsy] 
    % error_coef,: Maximum admitted value to validate the approximation of
    %            the imperfection. 
    % ModeType: Type of projected modes.
    %   'All': All modes
    % Nphi: Number of transverse modes
    % Lx, Ly: Plate dimensions
    
%% Output parameters:
    % proj: Vector containing the coefficients of the series expansion. 
    % modeIndices: kx and ky of the projected modes
    % Approximation: Profile obtained from the projection coefficients
    % error: Error of the approximation

    
%% %%%%%%%%%%%%%%%%%%%%%%%% Data obtention %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

[Kx, Ky] = meshgrid(1:Nphi/2, 1:Nphi/2);
Kx = reshape(Kx, [(Nphi/2)^2 1]);
Ky = reshape(Ky, [(Nphi/2)^2 1]);

Omega = ((Kx*pi/Lx).^2 + (Ky*pi/Ly).^2).^2;

[Omega, I] = sort(Omega);        

I(Nphi+1:end) = [];

modes = [Kx(I) Ky(I) Omega(I)];



switch ModeType
    case 'All'
        ValidModeIndices = 1:Nphi;
        
end

NumberOfModes = length(ValidModeIndices);

%% %%%%%%%%%%%%%%%%%%%% Initialization of variables %%%%%%%%%%%%%%%%%%%%%%

N=0;
 modeBase{Nphi}=[];
 modeIndices=zeros(Nphi,1);


[Nx, Ny] = size(Imperfection);

proj=zeros(Nphi,1); 

%%% Computation of the center of mass offset 
%zg=CartesianScalarProduct(Imperfection,ones(size(Imperfection)),0,Lx,0,Ly)/(Lx*Ly);
zg = 0;


Approximation=zg;

error=Inf;



%% %%%%%%%%%%%%%%%%%%%%%%%%% Main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while ((error>error_coef)&&(N<NumberOfModes))
       
    N=N+1;
    modeBase(N)=mat2cell(modes(ValidModeIndices(N),:),1); 
    
    modeIndices(N)=ValidModeIndices(N);
        
    %Projection of the imperfection
        
    [ ~, ~, phi] = ModeShapeRectangular( BC, modeBase{N}(1), modeBase{N}(2), Lx, Ly, Nx, Ny );

    phi = phi/sqrt(CartesianScalarProduct(phi,phi,0,Lx,0,Ly));

    proj(N) = CartesianScalarProduct(Imperfection,phi,0,Lx,0,Ly);
    
    Approximation = Approximation + proj(N)*phi;
    
    Approximation = Approximation - min(min(Approximation));
    
    %% Strict computation of the error. (Equivalent to the circular case)
    %auxImp = Imperfection;
    %auxImp(auxImp == 0) = 1;
    %error = max(std(abs(Imperfection-Approximation)./auxImp));
    
    %% Less strict computation of the error. 
    error = max(max(abs(Imperfection-Approximation)))/abs(max(max(Imperfection))-min(min(Imperfection)));        

end

proj=proj(1:N);

modeIndices=modeIndices(1:N);
 




