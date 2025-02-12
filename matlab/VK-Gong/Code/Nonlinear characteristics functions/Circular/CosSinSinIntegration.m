function [ S ] = CosSinSinIntegration( k, l, m )
% Computes the integration from 0 to 2pi of cos(k\theta)sin(l\theta)sin(m\theta)

% k, l and m can be matrices.

S = zeros(size(k));
I = find((k == abs(l-m) & l~= m & l~= 0 & m~= 0));
S(I) = pi/2*ones(size(I));
I = find(k == l+m & l~= 0 & m~= 0);
S(I) = -pi/2*ones(size(I));
I = find(k == 0 & l == m & l~= 0);
S(I) = pi*ones(size(I));

end

