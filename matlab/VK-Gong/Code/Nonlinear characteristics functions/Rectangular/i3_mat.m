function s = i3_mat( Npsi, Nphi, L)
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
            
            if n == 0 && p == 0
                
                s(m,n,p) = 0;
                
            elseif n == p && n ~= 0
                
                s(m,n,p) = -(- 4/L^3*(7*(-1)^(m1) + 8))*(L^4*(6*pi^2*p^2 - 2*pi^4*p^4))/(16*pi^4*p^4);
                
            else
                
                s(m,n,p) = (- 4/L^3*(7*(-1)^(m1) + 8))*(L*((6*L^3)/(n - p)^4 - (6*L^3)/(n + p)^4))/(2*pi^4) + (- 4/L^3*(7*(-1)^(m1) + 8))*(L*((3*L*cos(pi*(n + p))*(2*L^2 - L^2*pi^2*(n + p)^2))/(n + p)^4 - (3*L*cos(pi*(n - p))*(2*L^2 - L^2*pi^2*(n - p)^2))/(n - p)^4))/(2*pi^4);
                
            end
            
        end
        
    end
    
   
    
end
