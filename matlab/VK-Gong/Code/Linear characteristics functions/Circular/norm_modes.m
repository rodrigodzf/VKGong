function [Kkn] = norm_modes(k_t,xkn,R,dr, nu, KR, BC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%                                                                   %
%                              VK-Gong                              % 
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rr = 0:dr:R;
rr(1)=1e-10;

JJ0 = besselj(k_t,xkn*rr);
II0 = besseli(k_t,xkn*rr);


switch BC
    case {'free', 'elastic'}

        J0=besselj(k_t,xkn);
        J1=besselj(k_t-1,xkn);
        J2=besselj(k_t-2,xkn);

        I0=besseli(k_t,xkn); 
        I1=besseli(k_t-1,xkn);
        I2=besseli(k_t-2,xkn);

        Jtild=xkn^2.*J2+((nu-2*k_t+1)*xkn + KR)*J1+(k_t*(k_t+1)*(1-nu)-KR*k_t)*J0;
        Itild=xkn^2.*I2+((nu-2*k_t+1)*xkn + KR)*I1+(k_t*(k_t+1)*(1-nu)-KR*k_t)*I0;

        Rkn = JJ0 - (Jtild*II0/(Itild));
    
    case 'clamped'
        Jkn=besselj(k_t,xkn);
        Ikn=besseli(k_t,xkn);
        Rkn = JJ0*Ikn - Jkn*II0;
end

rR2 = rr.*Rkn.^2;

Kkn = sqrt(1./trapz(rr,rR2));

if k_t==0
    Kkn=Kkn/sqrt(2*pi);
else
    Kkn=Kkn/sqrt(pi);
end
%Rkn = Kkn*Rkn;





