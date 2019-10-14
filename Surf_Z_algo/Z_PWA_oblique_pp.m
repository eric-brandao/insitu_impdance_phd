function [R_PWA,alpha_PWA,Z_PWA]=Z_PWA_oblique_pp(k0,hs,z1,z2,r,p1,p2)

% R1=sqrt((r^2)+((hs-z1)^2));
% R1l=sqrt((r^2)+((hs+z1)^2));
% 
% R2=sqrt((r^2)+((hs-z2)^2));
% R2l=sqrt((r^2)+((hs+z2)^2));

R11=sqrt((r^2)+((hs-z1)^2));
R12=sqrt((r^2)+((hs+z1)^2));

R21=sqrt((r^2)+((hs-z2)^2));
R22=sqrt((r^2)+((hs+z2)^2));

Hw=p1./p2;

% R_PWA=(((exp(-1i*k0*R2)/R2)-Hw.*(exp(-1i*k0*R1)/R1))./...
%     ((exp(-1i*k0*R1l)/R1l).*Hw-(exp(-1i*k0*R2l)/R2l)));

% R_PWA=(((exp(-1i*k0*R21)/R21)-Hw.*(exp(-1i*k0*R11)/R11))./...
%     ((exp(-1i*k0*R12)/R12).*Hw-(exp(-1i*k0*R22)/R22)));

R_PWA=(((exp(-1i*k0*R11)/R11)-Hw.*(exp(-1i*k0*R21)/R21))./...
    ((exp(-1i*k0*R22)/R22).*Hw-(exp(-1i*k0*R12)/R12)));

alpha_PWA=1-((abs(R_PWA)).^2);

Z_PWA=(((1+R_PWA)./(1-R_PWA)).*(sqrt(r^2+hs^2)/hs).*...
    ((1i*k0*sqrt(r^2+hs^2))./((1i*k0*sqrt(r^2+hs^2))+1)));