function [R_PWA,alpha_PWA,Z_PWAc]=Z_PWA(k0,hs,z,r,Zm)

R1=sqrt((r^2)+((hs-z)^2));
R2=sqrt((r^2)+((hs+z)^2));

T1=((hs-z)/R1)*((1./(1i*k0*R1))+1);
T2=((hs+z)/R2)*((1./(1i*k0*R2))+1);
E=R2/R1;

R_PWA=(((Zm.*T1)-1)./((Zm.*T2)+1)).*E.*exp(-1i*k0*(R1-R2)); 
alpha_PWA=1-((abs(R_PWA)).^2);
Z_PWAc=(((1+R_PWA)./(1-R_PWA)).*(sqrt(r^2+hs^2)/hs).*...
    ((1i*k0*sqrt(r^2+hs^2))./((1i*k0*sqrt(r^2+hs^2))+1)));