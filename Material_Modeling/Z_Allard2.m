function Z=Z_Allard2(w,rho,dimm,resist,por,ks,s)

% Material data
% por=0.95; %porosidade
eta=1.84e-5; %viscosidade do ar
% a=5.41e-6; %raio das fibras
% ks=1; %1 para fibrosos
% s=1; % fator de estrutura
Np=0.77;   
gama=1.4;
Po=1.01320e5;

% Material geometry
%dimm=input('Espessura em mm: ');%6/1000; %espessura do material
% dimm=40;
di=dimm/1000;
%resist=input('Resistividade: ');%27.3*((1-por)^1.53)*(eta/(4*(a^2))) %eq 5.3
% resist=20000;
forma=(1/s)*sqrt((8*eta*ks)/(por*resist));

roe=((ks*rho).*(1+(((resist*por)./(1i.*w.*rho.*ks)).*...
    (sqrt(1+(((4*1i*(ks^2)*eta*rho).*w)./...
    ((resist^2)*(forma^2)*(por^2))))))));
Ke=((gama*Po)./(gama-((gama-1)./...
    (1+(((8*eta)./((1i*((2*forma)^2)*Np).*(w).*rho)).*...
    (sqrt(1+(((1i*rho).*(w).*(Np).*((2*forma)^2))./(16*eta)))))))));
zc=(sqrt(Ke.*roe));
kc=((w).*sqrt((roe)./(Ke)));
Z=-1i.*zc.*(cot(kc*di));
