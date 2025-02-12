function [ y ] = int1(m,p,L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Integral of Xdd(m,x)*Xdd(p,x) from zero to L, where X is the clamped plate
%function and d denotes derivative in x (dd = double derivative)

if m == 0 && p == 0
    
    
    y = 720/L^3;
    

elseif m == p 
    
                y = (pi^4*m^4 - 672*(-1)^m - 768)/(2*L^3);



elseif m == 0 || p == 0
    
    
        
                y = 0;


        
else

            y = -(24*(7*(-1)^m + 7*(-1)^p + 8*(-1)^m*(-1)^p + 8))/L^3;
            
end
        
    





