function s = i1_mat( Npsi, Nphi, L)
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
    m1 = m - 1;
    for n = 1 : Nphi 
        
        for p = 1 : Nphi 
            
            if n == 0 || p == 0
    
              s(m,n,p) = 0;
            
            elseif m1 == 0 && n == p && n ~= 0
                
                s(m,n,p) = L/2;
                
            elseif m1 == p - n || m1 == n - p
                
                s(m,n,p) = L/4;
                
            elseif m1 == - n - p || m1  == n + p
                
                s(m,n,p) = - L/4;
                
            end
            
        end
        
    end
    
end
            