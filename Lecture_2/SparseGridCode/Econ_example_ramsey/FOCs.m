%**************************************************************************
%
% This routine computes the firt-order conditions of the Ramsey model
%
%**************************************************************************

function [y]=FOCs(v,k,a,alpha,beta,delta,gamma,pers,wquad,pquad,z)

kp=real(v); 
ap=exp( (pers*log(a)*ones(1,size(pquad,1))+pquad)' );
kpvec=ones(size(ap,1),1)*kp;
statep = [kpvec ap];
kpp =spinterp(z,statep);
c=max(a*k^alpha + k*(1-delta) - kp,10^-6);
cp=max(kpvec.^alpha + kpvec*(1-delta)- kpp,10^-6);

% Euler error:
y = c / ( beta*wquad*((alpha*ap.*kpvec.^(alpha-1) + (1-delta)).*(cp).^(-gamma)))^(-1/gamma) - 1;
end

%**************************************************************************
