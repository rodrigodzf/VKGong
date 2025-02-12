function [ y ] = int4(m,p,L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integral of Xd(m,x)*Xd(p,x) from zero to L, where X is the clamped plate
%function and d denotes derivative in x

if m  == 0 && p == 0
    
    y = 120/(7*L);
    
        elseif m == p && m ~= 0
            
                y = (768*pi^2*m^2 - 47040*(-1)^m + 35*pi^4*m^4 + 432*(-1)^m*pi^2*m^2 - 53760)/(70*L*pi^2*m^2);
    
        elseif m == 0 
        
                y = (60*((-1)^p + 1)*(pi^2*p^2 - 42))/(7*L*pi^2*p^2);


        
        elseif p == 0
        
                y = (60*((-1)^m + 1)*(pi^2*m^2 - 42))/(7*L*pi^2*m^2);

        
else
    
                y = 192/35/L*(1 + (-1)^m*(-1)^p) - 192/m^2/p^2/L/pi^2*((p^2+m^2)*(1 + (-1)^m*(-1)^p)) -  168/m^2/p^2/L/pi^2*((p^2+m^2)*((-1)^m + (-1)^p)) + 108/35/L*((-1)^m + (-1)^p); 

end

end
        
    





