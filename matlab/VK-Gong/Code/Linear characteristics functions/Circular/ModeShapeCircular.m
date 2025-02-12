function [ U, V, W ] = ModeShapeCircular(k, c, xkn, R, nu, KR, kmax, nmax  )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Computes the transverse mode shapes of circular plate with boundary conditions corresponding to:
% Free edge: KR = 0
% Clamped edge: KR = inf;
% Elastic edge: KR

%% Inputs
% k: number of nodal radius
% n: number of nodal circles
% c: cos(1)/sin(2) mode 
% xkn: eigenfrequency
% KR: Normalized rotational stiffness (KT is considered when computing xkn)
% R: plate radius
% kmax: number of theta points
% nmax: number of points in x

%% Outputs
% U, V : Cartesian coordinates of points in W
% W : modal shape

    theta = linspace(0,2*pi, kmax);
    x = linspace(0, xkn, nmax);
    r=x*R/xkn;
    X=x'*ones(size(theta));
    THETA=ones(size(x))'*(theta);

      J=besselj(k,X);
      I=besseli(k, X);
      
      if KR == inf
          Jkn=besselj(k,xkn);
          Ikn=besseli(k,xkn);
          
          W=J-(Jkn*I/Ikn); 
          
      else
          J2=besselj(k-2,xkn);
          J1=besselj(k-1,xkn);
          J0=besselj(k,xkn);
          I2=besseli(k-2,xkn);
          I1=besseli(k-1,xkn);
          I0=besseli(k,xkn);
          Jtild=xkn^2.*J2+((nu-2*k+1)*xkn + KR)*J1+(k*(k+1)*(1-nu)-KR*k)*J0;
          Itild=xkn^2.*I2+((nu-2*k+1)*xkn + KR)*I1+(k*(k+1)*(1-nu)-KR*k)*I0;
      
          W=J - (Jtild*I/Itild); 
      end
       
      W=W.*cos(k*THETA+(c-1)/2*pi); 
    
      W=W/max(max(W)); 
      U=r'*cos(theta);
      V=r'*sin(theta);

end

