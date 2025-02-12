function s = i10_mat( Npsi, Nphi, L)
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
                
                s(m,n,p) = L^5/5;
                
            elseif n == p && p ~= 0
                
                s(m,n,p) = (15/L^4*((-1)^(m1) + 1))*(L^5*( + 4*pi^5*n^5 + 20*pi^3*n^3  - 30*pi*n))/(40*pi^5*n^5);
                
            else
                
                s(m,n,p) = -(15/L^4*((-1)^(m1) + 1))*(L*((4*pi*L^2*cos(pi*(n + p))*(6*L^2 - L^2*pi^2*(n + p)^2))/(n + p)^4 + (4*pi*L^2*cos(pi*(n - p))*(6*L^2 - L^2*pi^2*(n - p)^2))/(n - p)^4))/(2*pi^5);
                
            end
            
            
        end
        
    end
    
   
    
end