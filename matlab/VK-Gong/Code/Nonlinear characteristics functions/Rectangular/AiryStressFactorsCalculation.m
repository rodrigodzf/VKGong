function [coeff0, coeff1, coeff2] = AiryStressFactorsCalculation( BC, Npsi, Lx, Ly)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solves the eigenproblem associated to the Airy stress function for the
% inplane direction and returns the factors necessary for the computation
% of the H coupling coefficients. 

tic

switch BC            
    case {'SimplySupported'}

        str = sprintf('\n\n-------------------------eigen CALCULATION-------------------------\n\n');
        disp(str);
        
        %define indices (m,n) to use in loops

        [m,n] = meshgrid(0:Npsi-1, 0:Npsi-1);
        m = reshape(m, Npsi^2, 1);
        n = reshape(n, Npsi^2, 1);
        
        index_of_loop = [m n];


        %initialise stiffness and mass matrices
        K = zeros(Npsi^2,Npsi^2);
        M = zeros(Npsi^2,Npsi^2);






        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %             PART ONE, EIGENVALUE PROBLEM FOR AIRY STRESS FUNCTION

        for r = 1 : Npsi^2

            m = index_of_loop(r,1);

            n = index_of_loop(r,2);

            for s = 1 : Npsi^2

                p = index_of_loop(s,1);

                q = index_of_loop(s,2);

                K(r,s) =  int1(m,p,Lx)*int2(n,q,Ly) + int2(m,p,Lx)*int1(n,q,Ly) + 2*int4(m,p,Lx)*int4(n,q,Ly);
                M(r,s) =  int2(m,p,Lx)*int2(n,q,Ly);

            end

        end


        [VEC, VAL] = eig(K,M); %calculate eigenvalues and eigenvectors

        %find indices of negative or imaginary eigenvalues and remove them along
        %with corresponding eigenvectors
        [index_of_neg] = find(diag(VAL) < 0);
        [index_of_imag] = find(imag(diag(VAL)));
        VAL(:,index_of_neg) = [];
        VEC(:,index_of_neg) = [];
        VAL(:,index_of_imag) = [];
        VEC(:,index_of_imag) = [];

        %sort the eigenvalues and eigenvectors
        tres = find(VAL);
        [auto,v] = sort(VAL(tres));
        S = length(auto);


        coeff = VEC(:,v);

        coeff0 = zeros(Npsi^2,length(coeff(1,:)));
        coeff1 = zeros(Npsi^2,length(coeff(1,:)));
        coeff2 = zeros(Npsi^2,length(coeff(1,:)));

        % THIS ONE WORKS�����������������������������������������������������
        NN = int2_mat(Npsi,Lx);

        MM = int2_mat(Npsi,Ly);

        NN = reshape(NN,[Npsi^2,1]);

        MM = reshape(MM, [1, Npsi^2]);

        NN = repmat(NN,[1,Npsi^2]);

        MM = repmat(MM,[Npsi^2,1]);

        nmatr = NN.*MM;
        nmatr = full(nmatr);
        nmatr = reshape(nmatr,[Npsi Npsi Npsi Npsi]);
        nmatr = permute(nmatr,[4 1 3 2]);

        % nmatr = reshape(nmatr,[DIM*DIM*DIM DIM]);
        % nmatr = permute(nmatr,[2 1]);

        nmatr = reshape(nmatr,[1,Npsi^4]);
        nmatr = sparse(nmatr);

        %�����������������������������������������������������������������������




        for d=1:S

            temp = coeff(:,d);

            temp = reshape(temp, [Npsi^2,1]);

            temp = repmat(temp,[1 Npsi^2]);

            temp2 = permute(temp,[2,1]);

            temp3 = temp.*temp2;

            temp3 = reshape(temp3,[Npsi^4,1]);

            norms = (nmatr*temp3);

            coeff0(:,d) = coeff(:,d)/sqrt(norms);
            coeff1(:,d) = coeff(:,d)/sqrt(norms)/sqrt(auto(d));
            coeff2(:,d) = coeff(:,d)/sqrt(norms)/(auto(d));

        end


        S = floor(S/2);
        
        otherwise
            disp('Unknown boundary conditions');
end

coeff0 = coeff0(1:S,1:S);
coeff1 = coeff1(1:S,1:S);
coeff2 = coeff2(1:S,1:S);


toc