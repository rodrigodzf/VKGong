function [ fex ] = StrikeExcitation( T0, fm, Twid,  fs, Tn )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that generates a raised cosine time signal.

    fex = zeros(Tn,1);
    
    Nwid = round(Twid*fs);
    N0 = round(T0*fs);
    
    n = (1:(2*Nwid +1))';
        
    fex_temp = fm/2*(1 + cos(pi*((n-1-Nwid)/Nwid)));
    
    n = n+N0;
   
    if n(end)>Tn
        disp('Warning: Excitation lasts longer than the simulation time. Samples after Ts will be omitted.');
        disp(['Deleted samples: ', num2str(sum(n>Tn))]);
        n(n>Tn) = [];
        fex_temp = fex_temp(1:length(n));
    end
    
    
    fex(n) = fex_temp;




end

