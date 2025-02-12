function s = i4_mat( Npsi, Nphi, L)
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
                
                s(m+1,n+1,p+1) = 0;
                
                
            elseif n == p && n ~= 0
                
                s(m,n,p) = -(6/L^2*(2*(-1)^(m1) + 3))*(L^3*(6*pi*p - 4*pi^3*p^3))/(24*pi^3*p^3);
                
            else
                
                s(m,n,p) = (6/L^2*(2*(-1)^(m1) + 3))*(L^3*cos(pi*(n - p)))/(pi^2*(n - p)^2) - (6/L^2*(2*(-1)^(m1) + 3))*(L^3*cos(pi*(n + p)))/(pi^2*(n + p)^2);
                
            end
            
        end
        
        
    end
    
  
    
end