function [mode_l, Zeros  ] = ComputeInPlaneEigenfrequenciesCircular( dx, xmax, BC, nu )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function computes the inplane eigenvalues \zeta (\zeta^2 = \omega) 
% for a corresponding set of boundary conditions in a plate with unitary 
% radius. It also sorts the eigenvalues in terms of number of nodal radius 
% and circles. 


%% Input parameters
% dx: Precision step. 
% x_max: Top limit for \xi
% BC: Boundary conditions: "free", "clamped", "elastic". 
% nu: Poisson coefficient
% KR: Normalized rotational stiffness with respect to bending stiffness. KR = Kr/D 
% KT: Normalized transverse stiffness with respect to bending stiffness. KT = Kt/D;

%% Output parameters
% mode_l: Vector containing the transverse modes
%   <1:Index> <2:x_i> <3:k> <4:n> <5:c> <6:x^2_i>
% Zeros: Matrix containing the found zeros. 


switch BC
    case {'free', 'elastic'}
        R = 0;
        T = 0;
    case 'clamped'
        R = inf;
        T = inf;        
    otherwise
        disp('Unknown boundary conditions');
end


x=0:dx:0.0001*round(10000*xmax);

k=0; 

zer = 0;

Zeros = zer;

    while ~isempty(zer)
        switch BC            
            case {'free', 'elastic'}

                J0 = besselj(k-1,x);
                J1 = besselj(k,x);

                I0 = besseli(k-1,x);
                I1 = besseli(k,x);

                f = J1.*I0-J0.*I1;

            case 'clamped'
                J3=besselj(k-3,x);
                J2=besselj(k-2,x);
                J1=besselj(k-1,x);
                J0=besselj(k,x);

                I3=besseli(k-3,x);
                I2=besseli(k-2,x);
                I1=besseli(k-1,x);
                I0=besseli(k,x);

                Jtild=x.^2.*J2+(-nu-2*k+1)*x.*J1 + (k*(k+1) + nu*k*(1-k))*J0;
                Itild=x.^2.*I2+(-nu-2*k+1)*x.*I1 + (k*(k+1) + nu*k*(1-k))*I0;

                f=Itild.*(x.^3.*J3+(4-3*k)*x.^2.*J2+k*(k*(1-nu)-2)*x.*J1+(1+nu)*k^2*(1+k)*J0)-Jtild.*(x.^3.*I3+(4-3*k)*x.^2.*I2+k*(k*(1-nu)-2)*x.*I1+(1+nu)*k^2*(1+k)*I0);
        end
        
        zer = FindZeros(x,f);
        l = length(zer);
        Zeros(k+1,1:l) = zer;

        k = k+1;

    end
    
    mode_l = SortZeros(Zeros);
    
    k = 1;


    if ((R == 0) && (T == 0))
        mode_l(mode_l(:,3)~=1,4) =  mode_l(mode_l(:,3)~=1,4) -1;
    end


    
    switch BC
        case 'free'
            mode_l_free = mode_l;
            save('Parameters/Mode files/Circular/mode_l_free.mat', 'mode_l_free', 'Zeros');
        case 'clamped'
            mode_l_clamped = mode_l;
            save(sprintf('Parameters/Mode files/Circular/mode_l_clamped-nu_%d.mat', nu), 'mode_l_clamped', 'Zeros');
        case 'elastic'
            mode_l_elastic = mode_l;
            save('Parameters/Mode files/Circular/mode_l_elastic.mat', 'mode_l_elastic', 'Zeros');
        otherwise
            disp('Unknown boundary conditions');
    end
    
end

