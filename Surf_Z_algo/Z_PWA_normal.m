function [R_PWA,alpha_PWA,Z_PWA]=Z_PWA_normal(k0,hs,z,Zm)

T1=(1./(1i*k0*(hs-z)))+1;
T2=(1./(1i*k0*(hs+z)))+1;
E=(hs+z)/(hs-z);

R_PWA=(((Zm.*T1)-1)./((Zm.*T2)+1)).*E.*exp(2*1i*k0*z); 
alpha_PWA=1-((abs(R_PWA)).^2);
Z_PWA=(((1+R_PWA)./(1-R_PWA)).*((1i*k0*hs)./((1i*k0*hs)+1)));