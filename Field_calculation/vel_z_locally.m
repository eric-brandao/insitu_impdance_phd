function uz_loc=vel_z_locally(k0,hs,z,r,beta)
%%% This function calculates the sound field for a locally reactive
%%% infinite plane
%%% Uses adaptative quadrature integration
%%% Input parameters:
%%% k0 -wave number in air
%%% hs - source height
%%% z - sensor height
%%% r - horizontal distance
%%% beta - surface admitance: beta=1./Z
%%% References:
%%% [1] - X.DiandK.E.Gilbert,“An exact Laplace transform formulation 
%%%     for a point source above a ground surface," J.Acoust.Soc.Am.
%%%     93(2),pp.714-720,1993.

%% Geometry
R1=sqrt((r^2)+((hs-z)^2));
R2=sqrt((r^2)+((hs+z)^2));

%% Integrand
int_u=@(q)((exp((-q.*k0).*beta)).*...
    ((exp(-1i*k0.*(sqrt((r^2)+(hs+z-1i*q).^2))))./...
    (sqrt((r^2)+(hs+z-1i*q).^2))).*...
    ((hs+z-1i*q)./sqrt(r^2+(hs+z-1i*q).^2)).*...
    ((1./(1i*k0*(sqrt((r^2)+(hs+z-1i*q).^2))))+1));

%% Integral
Iq=quad(int_u,0,20);
uz_loc=((exp(-1i*k0*R1)/R1)).*((1./(1i*k0*R1))+1)*((hs-z)/R1)-...
    (exp(-1i*k0*R2)/R2).*((1./(1i*k0*R2))+1)*((hs+z)/R2)+...
    (2*k0*beta)*Iq;