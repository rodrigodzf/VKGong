function [c] = c_preset(X,om,Nphi, dFac, dExp, dCons)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function that loads the damping preset values. 
c = zeros(Nphi,1);



switch X 
    case 'Undamped'
        c = zeros(Nphi,1);
    case 'PowerLaw'
        for i = 1 : Nphi
            c(i) = dFac*om(i)^dExp + dCons;
        end
    otherwise 
        disp('Damping preset value not found.')
        
        
        
        
end

