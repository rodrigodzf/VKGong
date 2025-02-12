function [ y ] = int2_mat( N, L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%integral of X(m,x)*X(p,x) from zero to L, where X is the clamped plate
%function and d denotes derivative in x (dd = double derivative). The
%result is a matrix used to compute the norm of the eigenvectors. 

y = zeros(N,N);

for m = 0 : N - 1
    
    
    for p = 0 : N - 1
    

            y(m+1,p+1) = -(L*(11760*(-1)^m + 11760*(-1)^p - 16*pi^4*m^4 + 13440*(-1)^m*(-1)^p + (-1)^m*pi^4*m^4 + (-1)^p*pi^4*m^4 - 16*(-1)^p*pi^4*(m^4*(-1)^m) + 13440))/(70*pi^4*m^4) - (L*(13440*m^4 + 11760*(-1)^m*m^4 + 11760*(-1)^p*m^4 + 13440*(-1)^m*(-1)^p*m^4))/(70*pi^4*m^4*p^4);
    
    end
    y(m+1,m+1) = (67*L)/70 - ((-1).^m.*L)./35 - (768*L)./(pi^4*m.^4) - (672*(-1).^m.*L)./(pi^4*m.^4);
end

m = 0 : N - 1;


% 
 y(m+1,1) = (3*L*((-1).^m + 1).*(pi^4*m.^4 - 1680))./(14*pi^4*m.^4);
% 
 y(1,m+1) = (3*L*((-1).^m + 1).*(pi^4*m.^4 - 1680))./(14*pi^4*m.^4);
% 
 y(1,1) = (10*L)/7;
%     
% 
 end
%         
    





