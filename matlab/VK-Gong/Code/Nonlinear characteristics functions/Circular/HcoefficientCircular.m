function [ H ] = HcoefficientCircular( kp, kq, cp, cq, xip, xiq, ki, ci, zeta, nu, KR,  dr_H )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Computes H^i_pq tensor for a circular plate with BC boundary conditions.


%% Input parameters: 
% kp, kq, ki : Number of nodal diameters of p and q
% cp, cq, ci : Configuration of modes p and q. c=1 for cos, and c=2 for sin.
% xip, xiq : Squareroot of the angular eifrequency of modes p and q 
% zeta: Squareroot of the angular eifrequency of mode i
% nu: Poisson ratio
% KR: Normalized rotational stiffness KR = Kr/D, with Kr standing for the
% distributed rotational stiffness and D for the bending stiffness of the
% plate. 
% dr_H: Integration step for the computation of H. 

%% Output parameters
% H^i_{pq} = int_S \Psi_i L(\Phi_p,\Phi_q) dS





%% Conditions for (not) being null
% H = 0 if 
%     a) ki ~=  { kp + kq, abs(kp - kq) } or 
%     b) ki ==  { kp + kq, abs(kp - kq) } and the combination of c's is one
%     of the following:
%           cos / cos / sin
%           sin / sin / sin
%           sin / cos / cos
%           cos / sin / cos

I=[kp+kq abs(kp-kq)];
test=find(ki==I, 1);

if isempty(test)
  H=0;
  return
else
    if (((cp == cq) && (ci == 2)) || ((cp ~= cq) && (ci == 1)))
        H = 0;
        return
    end
end

rr=0:dr_H:1;
rr(1)=1e-8; %% To avoid singularities

%% ---------------------
% Calculation of R, R', R''
% ---------------------
% Transverse modes 

kk = [kp; kq];
xi = [xip; xiq];

R = zeros(2, length(rr));
dR = zeros(2, length(rr));
ddR = zeros(2, length(rr));

for ii = 1:2
  xkn = xi(ii);
  k = kk(ii);
  
  JJ0 = besselj(k,xkn*rr);
  JJ1 = besselj(k+1,xkn*rr);
  dJJ0 = -JJ1*xkn+k./rr.*JJ0;
  ddJJ0 = -(xkn^2+k./rr.^2-k^2./rr.^2).*JJ0+xkn./rr.*JJ1;
  
  II0 = besseli(k,xkn*rr);
  II1 = besseli(k+1,xkn*rr);
  dII0 = II1*xkn+k./rr.*II0;
  ddII0 = (xkn^2-k./rr.^2+k^2./rr.^2).*II0-xkn./rr.*II1;

  
  if KR==inf
      
      Jkn = besselj(k,xkn);
      Ikn = besseli(k,xkn);

      R(ii,:) = Ikn*JJ0-Jkn*II0;
      dR(ii,:) = Ikn*dJJ0-Jkn*dII0;
      ddR(ii,:) = Ikn*ddJJ0-Jkn*ddII0; 

  else
      J2 = besselj(k-2,xkn);
      J1 = besselj(k-1,xkn);
      J0 = besselj(k,xkn);
      I2 = besseli(k-2,xkn);
      I1 = besseli(k-1,xkn);
      I0 = besseli(k,xkn);
      Jtild = xkn^2.*J2+((nu-2*k+1)*xkn + KR)*J1+(k*(k+1)*(1-nu)-KR*k)*J0;
      Itild = xkn^2.*I2+((nu-2*k+1)*xkn + KR)*I1+(k*(k+1)*(1-nu)-KR*k)*I0;
    
      R(ii,:) = Itild*JJ0-Jtild*II0;
      dR(ii,:) = Itild*dJJ0-Jtild*dII0;
      ddR(ii,:) = Itild*ddJJ0-Jtild*ddII0;
      
  end
       
  % --------------------------------------   
  % Normalization of R (int(phi^2 dS)=1)
  % --------------------------------------

  
  rR2 = rr.*R(ii,:).^2;
  
  Kkn = sqrt(1./trapz(rr,rR2)); %% Multiplied by rr for the polar coordinates integration and squared to compute the norm.
  
  if k==0
    Kkn = Kkn/sqrt(2*pi);
  else
    Kkn = Kkn/sqrt(pi);
  end
  R(ii,:) = Kkn*R(ii,:);
  dR(ii,:) = Kkn*dR(ii,:);
  ddR(ii,:) = Kkn*ddR(ii,:);

end

%% ---------------------
% Calculation of S, S', S''
% ---------------------
% In-plane modes 

  if KR==inf
    J2  = besselj(ki-2,zeta);
    J1 = besselj(ki-1,zeta);
    J0 = besselj(ki,zeta);

    I2 = besseli(ki-2,zeta);
    I1 = besseli(ki-1,zeta);
    I0 = besseli(ki,zeta);

    Jtild = zeta.^2.*J2+(-nu-2*ki+1)*zeta.*J1+(ki*(ki+1)+nu*ki*(1-ki))*J0;
    Itild = zeta.^2.*I2+(-nu-2*ki+1)*zeta.*I1+(ki*(ki+1)+nu*ki*(1-ki))*I0;   

    J = besselj(ki,zeta*rr);
    I = besseli(ki,zeta*rr);

    S = Itild*J-Jtild*I;


  else

    J =  besselj(ki,zeta);
    I = besseli(ki,zeta);

    JJ0 = besselj(ki,zeta*rr);


    II0 = besseli(ki,zeta*rr);


    S = I.*JJ0-J.*II0;

  end
  
    % --------------------------------------
    % normalization of S (int(psi^2 dS)=1)
    % --------------------------------------

    rS2 = rr.*S.^2; %% Multiplied by rr for the polar coordinates integration and squared to compute the norm.
    Lkn = sqrt(1./trapz(rr,rS2));

    if ki==0
        Lkn = Lkn/sqrt(2*pi);
    else
        Lkn = Lkn/sqrt(pi);
    end
    
    S = Lkn.*S;

    
    %% Computation of the coefficients

    fctH1 = S.*(ddR(1,:).*(dR(2,:)-kq^2*R(2,:)./rr)+ddR(2,:).*(dR(1,:)-kp^2*R(1,:)./rr));
    fctH2 = S.*(dR(1,:)-R(1,:)./rr).*(dR(2,:)-R(2,:)./rr)./rr;

    % Radius dependent term
    H1 = trapz(rr,fctH1);
    H2 = trapz(rr,fctH2);

    % Theta dependent term
    if cp==1
      if cq==1
        beta1 = CosCosCosIntegration(kp,kq,ki);
        beta2 = CosSinSinIntegration(ki,kp,kq);
      else
        beta1 = CosSinSinIntegration(kp,kq,ki);
        beta2 = -CosSinSinIntegration(kq,kp,ki);
      end
    else
      if cq==1
        beta1 = CosSinSinIntegration(kq,kp,ki);
        beta2 = -CosSinSinIntegration(kp,kq,ki);
      else
        beta1 = CosSinSinIntegration(ki,kp,kq);
        beta2 = CosCosCosIntegration(kq,kp,ki);
      end
    end
    
    
    H = H1.*beta1-2*kp*kq*H2.*beta2;

end

