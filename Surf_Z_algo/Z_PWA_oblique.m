function [R_PWA,alpha_PWA,Z_PWA]=Z_PWA_oblique(k0,hs,z,r,Zm)

%%%%%%%%%%%%%%%% Input list %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k0 --> wave number vector; vector of the same length of f
% hs --> Source height: vertical distance from the monopole source to the 
%        surface of the sample
% z --> Sensor height: same definition as hs
% r --> Source to sensor horizontal separation (if they are on the same
%       line r=0, which means normal incidence)
% Zm --> The impedance measured by the PU sensor (e.g 1 cm above the
%        surface of the sample); vector of the same length
%        of f

%%%%%%%%%%%%%%%% Output list %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R_PWA --> Reflection coefficient calculated with the image source method
% alpha_PWA --> Absoption coefficient calculated with the image source method
% Z_PWA --> Surface impedance calculated with the image source
%           method; first estimative to q-term iterative algorithm

R1=sqrt((r^2)+((hs-z)^2));
R2=sqrt((r^2)+((hs+z)^2));

T1=((hs-z)/R1)*((1./(1i*k0*R1))+1);
T2=((hs+z)/R2)*((1./(1i*k0*R2))+1);
E=R2/R1;

R_PWA=(((Zm.*T1)-1)./((Zm.*T2)+1)).*E.*exp(-1i*k0*(R1-R2)); 
alpha_PWA=1-((abs(R_PWA)).^2);
Z_PWA=(((1+R_PWA)./(1-R_PWA)).*(sqrt(r^2+hs^2)/hs).*...
    ((1i*k0*sqrt(r^2+hs^2))./((1i*k0*sqrt(r^2+hs^2))+1)));