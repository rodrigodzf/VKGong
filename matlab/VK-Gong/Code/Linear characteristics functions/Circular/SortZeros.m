function [ TAB ] = SortZeros( Zeros)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This function sorts the zeros in matrix Zeros in ascending order. 

%% Input parameters
% Zeros contains the zeros obtained gradually. (k row, n column)

%% Output parameters
% TAB: Sorts the zeros in a matrix with the following column order
% <1:Index> <2:x_i> <3:k> <4:n> <5:c> <6:x^2_i>

TAB=[];
ii=1;

Zeros(Zeros==0)=Inf*ones(size(find(Zeros==0)));
zer=min(min(Zeros));

while zer~=Inf
  [k,n]=find(Zeros==zer);
  TAB(ii,:)=[ii zer k-1 n 1 zer^2 ];        
  ii=ii+1;  
  if k>1
      TAB(ii,:)=[ii zer k-1 n 2 zer^2 ];
      ii=ii+1;
  end
  Zeros(k,n)=Inf;
  zer=min(min(Zeros));    
end




end

