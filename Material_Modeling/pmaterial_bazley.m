function [Z Zc kc]=pmaterial_bazley(c0,rho0,f,resist,dimm)
%%% This function calculates the characteristic impedance, complex wave
%%% number, complex density and the surface impedance for a porous sample
%%% of thickness d1 [mm] backed by an impervious rigid wall.

%%% Reference
%%% [1] - J.F. Allard and N. Atalla, "Propagation of sound in porous media:
%%%            modeling sound absorbing materials",Chapter 2,pp.20-23 Wiley, 2009.

%%% [2] - T.J. Cox, p. D’Antonio and M. Schroeder, "Acoustic Absorbers and
%%% Diffusers, Theory, design and application", Chapter 5, pp. 173-174,
%%% ASA, 2005

%% Input
di=dimm/1000;                                   % Thickness of the sample [m]
w=2*pi*f; k0=w/c0;
%% Calculations
X=rho0*f/resist; 
minX=rho0*f(1)/resist;
maxX=rho0*f(length(f))/resist;
% Zc=(rho0*c0)*(1+(0.0571*(X.^-0.754))-(1i*0.087*(X.^-0.732)))
% kc=k0.*(1+(0.0978*(X.^-0.700))-(1i*0.189*(X.^-0.595)));

X=1000*f/resist;
% Zc=((rho0*c0)*(1+(9.08*(X.^-0.754))-(1i*11.9*(X.^-0.732))));
Zc=((rho0*c0)*(1+(9.08*(X.^-0.75))-(1i*11.9*(X.^-0.73))));
% kc=-1i*(k0.*(1+(10.8*(X.^-0.700))-(1i*10.3*(X.^-0.595))));
% kc=1i*(k0.*(10.3*(X.^-0.595)+1i*(1+(10.8*(X.^-0.700)))));
kc=-1i*(k0.*(10.3*(X.^-0.59)+1i*(1+(10.8*(X.^-0.700)))));

Z=-1i.*Zc.*(cot(kc*di));                        % Surface impedance