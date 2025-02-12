function s = i5_mat( Npsi, Nphi, L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary integral for the computation of the coupling coefficient H

s = zeros(Npsi,Nphi,Nphi);

for m = 1 : Npsi 
    %m1 = m - 1;
    for n = 1 : Nphi 
        
        for p = 1 : Nphi 
            
            if n == p && n ~= 0
                
                s(m,n,p) = L/2;
                
                
            end
            
        end
        
    end
    
end
    
    s = -s;