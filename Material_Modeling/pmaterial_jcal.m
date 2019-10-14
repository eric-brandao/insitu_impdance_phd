function [Z, Zc, kc, rho_e]=pmaterial_jcal(w,rho0,dimm,resist,por,a_inf,Lam,Lam_l)
%%% This function calculates the characteristic impedance, complex wave
%%% number, complex density and the surface impedance for a porous sample
%%% of thickness d1 [mm] backed by an impervious rigid wall.

%%% inputs:
% w - anfular frequency
% rho0 - air density
% dimm - thickness in mm
% resist - Flow resistivity
% por - porosity
% a_inf - tortuosity
% Lam - Vicous characteristic length
% Lam_l - Thermal characteristic length

%%% outputs
% Z - surface impedance
% Zc - characteristic impedance
% kc - complex wave number
% rho_e - characteristic density

%%% Reference
%%% [1] - J.F. Allard and N. Atalla, "Propagation of sound in porous media:
%%%            modeling sound absorbing materials", Wiley, 2009.

%% Constant Parameters
eta=1.84e-5;                    % Air viscosity
B2=0.77;                        % Prandtl Number   
gama=1.4;                       % cp/cv ratio
P0=1.01320e5;                   % Startic pressure
v=eta/rho0;
v_l=v/B2;

%% Input parameters
di=dimm/1000;                  % Thickness of the sample [m]

%% Calculations
q0=eta./resist;
q0_l=por.*(Lam_l.^2)/8;

Gw=sqrt(1+(((2*a_inf.*q0)./(por.*Lam)).^2).*(1i*w/v));
% Gw=sqrt(1+(((2*a_inf*q0)./(por*Lam))).*(1i*w/v));
Gw_l=sqrt(1+((Lam_l/4).^2).*(1i*w./v_l));

rho_e=(rho0*(a_inf+((v*por)./(1i*w*q0)).*Gw));  %Effective density

K_e=(gama*P0./(gama-((gama-1)./...
    (1+(((v_l*por)./(1i*w*q0_l)).*Gw_l)))));    % Effective bulk modullus

Zc=sqrt(K_e.*rho_e);                            % Characteristic impedance
kc=(w).*sqrt((rho_e)./(K_e));                   % Wave number
Z=-1i.*Zc.*(cot(kc.*di));                        % Surface impedance
