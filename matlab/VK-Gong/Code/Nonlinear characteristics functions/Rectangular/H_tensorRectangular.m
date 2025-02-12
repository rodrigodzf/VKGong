function [ H0, H1, H2 ] = H_tensorRectangular( coeff0, coeff1, coeff2, Nphi, Npsi, Lx, Ly, mode_t, BC )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of the H matrix

tic


str = sprintf('\n\n-------------------------H CALCULATION-------------------------\n\n');
disp(str);


S = length(coeff1(1,:));

% Initialise H tensor
H0 = zeros(S,Nphi*Nphi);
H1 = zeros(S,Nphi*Nphi);
H2 = zeros(S,Nphi*Nphi);

switch BC
    case 'SimplySupported'
        %call functions for H calculation (integral definitions!)
        m1 = g1(Npsi,Nphi,S,Lx, mode_t);
        m2 = g2(Npsi,Nphi,S,Lx, mode_t);
        m3 = g3(Npsi,Nphi,S,Ly, mode_t);
        m4 = g4(Npsi,Nphi,S,Ly, mode_t);
        m5 = g5(Npsi,Nphi,S,Lx, mode_t);
        m6 = g6(Npsi,Nphi,S,Ly, mode_t);

        for n = 1 : S

             f0 = coeff0(:,n).'*(m1.*m4 + m2.*m3 - 2*m5.*m6);
             f1 = coeff1(:,n).'*(m1.*m4 + m2.*m3 - 2*m5.*m6);
             f2 = coeff2(:,n).'*(m1.*m4 + m2.*m3 - 2*m5.*m6);

             H0(n,:) = f0;
             H1(n,:) = f1;
             H2(n,:) = f2;

        end

        % Put constants in
        H0 = H0*4*pi^4/Lx^3/Ly^3;
        H1 = H1*4*pi^4/Lx^3/Ly^3;
        H2 = H2*4*pi^4/Lx^3/Ly^3;

        for ii = 1 : S

            temp0 = H0(ii,:);
            temp1 = H1(ii,:);
            temp2 = H2(ii,:);

            v0 = find(abs(temp0/max(abs(temp0))) < 1e-8);
            v1 = find(abs(temp1/max(abs(temp1))) < 1e-10);
            v2 = find(abs(temp2/max(abs(temp2))) < 1e-8);

            temp0(v0) = 0; 
            temp1(v1) = 0; 
            temp2(v2) = 0; 

            H0(ii,:) = temp0;
            H1(ii,:) = temp1;
            H2(ii,:) = temp2;

        end
        
        otherwise
            disp('Unknown boundary conditions');
end

toc



end

