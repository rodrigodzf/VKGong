function s = i13_mat( Npsi, Nphi, L)
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
   
    for n = 1 : Nphi 
        
        for p = 1 : Nphi             
            if n == 0 && p == 0
                
                s(m,n,p) = L;
                
            elseif n == p
                
                s(m,n,p) = L/2;
                
            end
            
        end
        
    end
    
end


s = -s;