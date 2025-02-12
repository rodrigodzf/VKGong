function [ m ] = g5( Npsi, Nphi, S, Lx, mode_t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integration term included in the computation of the H coefficients.

s = zeros(Npsi,Nphi,Nphi);
s1 = i9_mat(Npsi,Nphi,Lx);
s2 = i10_mat(Npsi,Nphi,Lx);
s3 = i11_mat(Npsi,Nphi,Lx);
s4 = i12_mat(Npsi,Nphi,Lx);
s5 = i13_mat(Npsi,Nphi,Lx);



for m = 1 : Npsi
    
    for n = 1 : Nphi
        
        n1 = mode_t(n,2); 
        
        for p = 1 : Nphi
            
             p1 = mode_t(p,2);
             
             s(m,n,p) = s1(m,n1,p1) + s2(m,n1,p1) + s3(m,n1,p1) + s4(m,n1,p1) + s5(m,n1,p1);
             s(m,n,p) = s(m,n,p)*p1*n1;
            
        end
        
    end
    
end



s = reshape(s,[Npsi Nphi^2]);


m = zeros(Npsi^2,Nphi^2);

for u = 1 : Nphi^2
    for i = 1 : Npsi
        m(Npsi*(i-1) + 1 :  Npsi*i, u) = s(i,u);
    end
end
m = m(1:S,:);
