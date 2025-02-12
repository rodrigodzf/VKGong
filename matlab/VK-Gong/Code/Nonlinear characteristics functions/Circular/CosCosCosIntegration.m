function [ S ] = CosCosCosIntegration( k, l, m )

% Computes the integration from 0 to 2pi of cos(k\theta)cos(l\theta)cos(m\theta)

% k, l and m can be matrices.

S=zeros(size(k));
I=find(m==l+k | m==abs(l-k)); % | l == m+k | l == abs(m-k) | k == l+m | k == abs(l- m));
S(I) = pi/2*ones(size(I));
I = find((k == 0 & l == m) | (l == 0 & k == m) |(m == 0 & k == l));
S(I) = pi*ones(size(I));
I = find(k == 0 & l == 0 & m == 0);
S(I) = 2*pi*ones(size(I));


end

