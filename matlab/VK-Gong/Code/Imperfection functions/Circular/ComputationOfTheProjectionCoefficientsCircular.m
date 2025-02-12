function [proj,modeIndices,Approximation, error] = ComputationOfTheProjectionCoefficientsCircular( Imperfection,error_coef,nu, KR, ModeType, mode_t  )
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
%	Imperfection: Values of the imperfection profile in points $[U,V]$.
%   error_coef: Top error admitted in the approximation of the imperfection. $error\_coef \in \left[0,1\right]$ 
%   nu: Poisson ratio. 
%   KR: Normalized rotational stiffness. Only used when BC = 'elastic'.
%   ModeType: Type of modes considered in the approximation of the imperfection. Possible values: 'All', 'Axisymmetric'.
%   mode_t: Vector containing the information corresponding to the translational vibration modes. For every mode, \\	
    
%% Output parameters:
    % proj: Vector containing the coefficients of the series expansion. 
    % Approximation: Imperfection approximated by the projection coefficients    
    % modeIndices: Indices of modes used for the approximation
    % error: Error in the approximation


    
%% %%%%%%%%%%%%%%%%%%%%%%%% Data obtention %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    


switch ModeType
    case 'All' %% All modes are considered. 
        ValidModeIndices=(1:size(mode_t,1));
        
    case 'Axisymmetric' %% Only the axysimmetric modes are considered, i.e. k=0.
        ValidModeIndices=find(mode_t(:,3)==0);

end
        
        
NumberOfModes=length(ValidModeIndices);

%% %%%%%%%%%%%%%%%%%%%% Initialization of variables %%%%%%%%%%%%%%%%%%%%%%
a=1; % Normalized radius used to compute the modeshapes. 

N=0;
modeBase{NumberOfModes}=[];
modeIndices=zeros(NumberOfModes,1);


[r,theta]=size(Imperfection);
kmax=theta; 
nmax=r;

proj=zeros(NumberOfModes,1); 

%%% Computation of the center of mass offset 
zg=PolarScalarProduct(Imperfection,ones(size(Imperfection)),0,1,0,2*pi)/(pi*a^2); %% This is 0 for the spherical cap, and should be checked for other axysimmetric imperfections. But it should be zero as well since y0 is included in the computation. 

Approximation=zg;

error=Inf;



%% %%%%%%%%%%%%%%%%%%%%%%%%% Main loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while ((error>error_coef)&&(N<NumberOfModes))
       
    N=N+1;
    modeBase(N)=mat2cell(mode_t(ValidModeIndices(N),2:5),1); 
    
    modeIndices(N)=ValidModeIndices(N);
        
    %Projection of the imperfection
            xkn = modeBase{N}(1);
            k = modeBase{N}(2);
            c=modeBase{N}(4);

            [~,~,phi] = ModeShapeCircular(k, c, xkn, a, nu, KR, kmax, nmax  );
            phi=phi/sqrt(PolarScalarProduct(phi,phi,0,1,0,2*pi));

            proj(N)=PolarScalarProduct(Imperfection,phi,0,1,0,2*pi);
            Approximation=Approximation+proj(N)*phi;

        error=max(std(abs(Imperfection-Approximation)./Imperfection));
            

end


proj=proj(1:N);

modeIndices=modeIndices(1:N);
 




