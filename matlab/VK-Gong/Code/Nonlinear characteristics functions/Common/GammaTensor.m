function [ G ] = GammaTensor( H1, filename, Nphi, Npsi )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function computes the Gamma tensor from H1 = H/\zeta_i.

% Input parameters
% H1 = H^i_pq/\zeta^2_i
% filename: Name of the output file. 
% Nphi: Number of transverse modes
% Npsi: Number of inplane modes

% Output parameters
% G: Gamma tensor size(G) = [Nphi Nphi Nphi Nphi]

tic

%Transform H tensor into a matrix and impose sparse description
H1 = sqrt(0.5)*H1(1:Npsi,1:Nphi,1:Nphi); %% The factor sqrt(0.5) is added to fulfill the Gamma definiton.
H1 = reshape(H1,[Npsi Nphi*Nphi]);

%Calculate Gamma in term of H
G = H1'*H1;

%Permute dimensions to get the coefficients in the right place
G = reshape(G,[Nphi Nphi Nphi Nphi]);
G = permute(G, [4 1 2 3]);
G = reshape(G,[Nphi*Nphi Nphi*Nphi]); 

G = full(G);
G = reshape(G,[Nphi Nphi Nphi Nphi]);

toc

save(filename, 'G');


end

