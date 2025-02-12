function [ m ] = g4( Npsi, Nphi, S, Ly, mode_t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integration term included in the computation of the H coefficients.

s = zeros(Npsi,Nphi,Nphi);
s1 = i1_mat(Npsi,Nphi,Ly);
s2 = i2_mat(Npsi,Nphi,Ly);
s3 = i3_mat(Npsi,Nphi,Ly);
s4 = i4_mat(Npsi,Nphi,Ly);
s5 = i5_mat(Npsi,Nphi,Ly);



for m = 1 : Npsi
    
    for n = 1 : Nphi
        
        n1 = mode_t(n,3); 
        
        for p = 1 : Nphi
            
             p1 = mode_t(p,3);
             
             s(m,n,p) = s1(m,n1,p1) + s2(m,n1,p1) + s3(m,n1,p1) + s4(m,n1,p1) + s5(m,n1,p1);
             s(m,n,p) = s(m,n,p)*(p1)^2;
            
        end
        
    end
    
end



s = reshape(s,[Npsi Nphi^2]);


m = repmat(s,[Npsi,1]);
m = m(1:S,:);
