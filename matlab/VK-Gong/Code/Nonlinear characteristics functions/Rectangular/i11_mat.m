function s = i11_mat( Npsi, Nphi, L)
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
                
                s(m,n,p) = (- 4/L^3*(7*(-1)^(m1) + 8))*L^4/4;
                
            elseif n == p && p ~= 0
                
                s(m,n,p) = (- 4/L^3*(7*(-1)^(m1) + 8))*L^4/8 + (- 4/L^3*(7*(-1)^(m1) + 8))*(3*L^4)/(8*pi^2*p^2);
                
            else
                
                s(m,n,p) = (- 4/L^3*(7*(-1)^(m1) + 8))*(L*((6*L^3)/(n - p)^4 + (6*L^3)/(n + p)^4))/(2*pi^4) - (- 4/L^3*(7*(-1)^(m1) + 8))*(L*((3*L*cos(pi*(n + p))*(2*L^2 - L^2*pi^2*(n + p)^2))/(n + p)^4 + (3*L*cos(pi*(n - p))*(2*L^2 - L^2*pi^2*(n - p)^2))/(n - p)^4))/(2*pi^4);
                
            end
            
        end
        
    end
    
end

