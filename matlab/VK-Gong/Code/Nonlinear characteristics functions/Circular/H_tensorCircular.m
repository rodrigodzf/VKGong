function [ H0, H1, H2 ] = H_tensorCircular( Nphi, Npsi, BC, nu, KR, KT, dr_H )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This functions compute the H matrices for a circular plate of radius 
% equal to unity with the boundary conditions stated by BC. 

%% Input parameters: 
% Nphi: Number of transverse modes
% Npsi: Number of in-plane modes
% BC: Boundary conditions at the edge: 'free', 'clamped', 'elastic'
% nu: Poisson coefficient
% KR: Normalized rotational stiffness KR = Kr/D, with Kr standing for the
% distributed rotational stiffness and D for the bending stiffness of the
% plate. 
% KT: Normalized translational stiffness KT = Kt/D, with Kt standing for
% the distributed rational stiffness at the edge. 
% dr_H: Integration step for the computation of H. 

%% Output parameters
% H0 = H^i_pq
% H1 = H0/zeta_i^2
% H2 = H0/zeta_i^4

% H(i,p,q)
%
% i : Index of in-plane mode i
% p : Index of transverse mode p
% q : Index of transverse mode q


switch BC
    case 'free'
        KR = 0;
        modeTfilename = sprintf('mode_t_free-nu_%f.mat', nu);
        load(modeTfilename);
        load mode_l_free
        mode_t = mode_t_free;
        clear mode_t_free
        mode_l = mode_l_free;
        clear mode_l_free;
        filename = sprintf('Parameters/H files/Circular/H_free-Nphi_%d-Npsi_%d-nu_%f-dr_%f.mat',Nphi,Npsi,nu,dr_H);
        
        
    case 'clamped'
        KR = inf;
        load mode_t_clamped
        modeLfilename = sprintf('mode_l_clamped-nu_%d.mat', nu);
        load(modeLfilename);
        mode_t = mode_t_clamped;
        clear mode_t_clamped
        mode_l = mode_l_clamped;
        clear mode_l_clamped
        filename = sprintf('Parameters/H files/Circular/H_clamped-Nphi_%d-Npsi_%d-nu_%f-dr_%f.mat',Nphi,Npsi,nu,dr_H);
        
    case 'elastic'
        modeTfilename = sprintf('mode_t_elastic-nu_%f-KR_%f-KT_%f.mat', nu, KR, KT);
        load(modeTfilename);
        load mode_l_free
        mode_t = mode_t_elastic;
        clear mode_t_elastic
        mode_l = mode_l_free;
        clear mode_l_free;
        filename = sprintf('Parameters/H files/Circular/H_elastic-Nphi_%d-Npsi_%d-nu_%f-dr_%f-KR_%f-KT_%f.mat',Nphi,Npsi,nu,dr_H,KR,KT);
        
    otherwise
        disp('Unknown boundary conditions');
end

% Transverse modes
k_t = mode_t(1:Nphi,3); % Number of nodal diameters
%n_t = mode_t(1:Nmax,4);
c_t = mode_t(1:Nphi,5); % cos(1)/sin(2) mode
xi_t = mode_t(1:Nphi, 2); 

% In-plane modes
k_l = mode_l(1:Npsi,3); % Number of nodal diameters
%n_l = mode_l(1:Npsi,4);
c_l = mode_l(1:Npsi,5); % cos(1)/sin(2) mode
zeta_l = mode_l(1:Npsi,2);

H0 = zeros(Npsi, Nphi,Nphi);
H1 = zeros(Npsi, Nphi,Nphi);
H2 = zeros(Npsi, Nphi,Nphi);

tic

for p = 1:Nphi
    for q = p:Nphi
        for i = 1:Npsi
            
            H0(i,p,q) = HcoefficientCircular( k_t(p), k_t(q), c_t(p), c_t(q), xi_t(p), xi_t(q), k_l(i), c_l(i), zeta_l(i), nu, KR, dr_H );
            
            H0(i,q,p) = H0(i,p,q);
            
            H1(i,p,q) = H0(i,p,q)/zeta_l(i)^2;
            
            H1(i,q,p) = H0(i,q,p)/zeta_l(i)^2;
            
            H2(i,p,q) = H0(i,p,q)/zeta_l(i)^4;
            
            H2(i,q,p) = H0(i,q,p)/zeta_l(i)^4;
        
        end      
        
    end
end

toc


save(filename, 'H0', 'H1', 'H2', '-v7.3');
end

