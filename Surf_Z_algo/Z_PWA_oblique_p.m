function [R_PWA,alpha_PWA,Z_PWA]=Z_PWA_oblique_p(k0,hs,z,r,Vp)

R1=sqrt((r^2)+((hs-z)^2));

R2=sqrt((r^2)+((hs-z)^2));

R_PWA=Vp*(R1/R2).*(exp(-1i*k0*(R2-R1)));

alpha_PWA=1-((abs(R_PWA)).^2);

Z_PWA=(((1+R_PWA)./(1-R_PWA)).*(sqrt(r^2+hs^2)/hs).*...
    ((1i*k0*sqrt(r^2+hs^2))./((1i*k0*sqrt(r^2+hs^2))+1)));