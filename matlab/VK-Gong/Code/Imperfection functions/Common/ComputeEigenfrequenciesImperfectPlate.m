function [Omega , A_matrix] = ComputeEigenfrequenciesImperfectPlate(proj,NA ,modeIndices, om, e, GammaFileName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
Nproj=length(proj);

%% Load Gamma coefficients Gamma^u_rps
load(GammaFileName);

gamma=G(1:NA,1:NA,1:NA,1:NA);
om = om(1:NA);

clear G;

ai=zeros(NA,1);

ai(modeIndices)=proj;

if length(ai)>NA
    disp('Eigenfrequencies computation: Number of projected modes is higher than computed eigenfrequencies');
    ai=ai(1:NA);
    modeIndices(modeIndices>NA)=[];
    Nproj=length(modeIndices);
end


air=zeros(NA,NA,NA,NA);
ais=zeros(NA,NA,NA,NA);
for i=1:Nproj
    air(:,modeIndices(i),:,:)=ai(modeIndices(i));
    ais(:,:,:,modeIndices(i))=ai(modeIndices(i));
end

alpha=2*e*gamma.*air.*ais;

clear gamma;

alpha=sum(alpha,4);
alpha=sum(alpha,2);
alpha=reshape(alpha,NA,NA);

A_matrix=alpha+diag(om.^2);
Omega=sort(sqrt(eig(A_matrix)));



    
end

