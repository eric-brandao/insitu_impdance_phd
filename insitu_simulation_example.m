clc; clearvars; close all
%%%% A minimalistic example on the simulation of in situ impedance
%%%% measurement.

%% Set some general configurations
general.freq = 100:10:10000; %%% frequency vector
general.c0 = 343; %%% sound speed
general.rho0 = 1.21; %%% air density

%% Porous material simulation
porous.resist = 9209;
porous.porosity = 0.99;
porous.a_inf = 1.00;
porous.Lam = 300e-6;
porous.Lam_l = 600e-6;
porous.thickness = 40/1000;

[porous.Z, ~, ~, ~]=pmaterial_jcal(2*pi*general.freq,general.rho0,...
    porous.thickness*1000,porous.resist,porous.porosity,...
    porous.a_inf,porous.Lam,porous.Lam_l);
porous.Vp = (porous.Z-general.c0*general.rho0)./...
    (porous.Z+general.c0*general.rho0);
porous.alpha = 1 - (abs(porous.Vp)).^2;

%%
%%% let us simulate the acoustic field. We'll define a sound source, 
%%% receivers and we will calculate the sound field for posterior
%%% processing. This emulates what the sensors will capture.
%%% source
source.hs = 0.3; %% source height
source.r = 0.0; %% horizontal distance between source-receiver
%%% Receivers
receiver.hr = 0.01;
%%% sound field simulation
measurement(1).p = zeros(1, length(general.freq));
measurement(1).uz = zeros(1, length(general.freq));
hq = waitbar(0, 'Calculating Field...');
for jf = 1:length(general.freq)
    waitbar(jf/length(general.freq),hq)
    k0 = 2 * pi * general.freq(jf)/general.c0;
    beta = (general.rho0*general.c0) / porous.Z(jf);
    measurement.p(jf) = pres_locally(k0,source.hs,receiver.hr,source.r,beta);
    measurement.uz(jf) = vel_z_locally(k0,source.hs,receiver.hr,source.r,beta);
end
close(hq)
clear jf beta k0 hq
%%%% The impedance measured by PU sensor
measurement.Zm = measurement.p./measurement.uz;

%% Now, let us recover the surface impedance.
%%%% First, assuming a specular reflection of spherical wave
[~,alpha_PWA, Z_PWA]=Z_PWA_oblique(2*pi*general.freq/general.c0,...
    source.hs,receiver.hr,source.r,measurement.Zm);

%%% Now the q-term model
Zq = Z_qterm_quad_2(general.freq,general.c0,source.hs,...
    receiver.hr,source.r,Z_PWA , measurement.Zm);
alpha_q = 1 - (abs((Zq - 1)./(Zq + 1))).^2;
%% plots
figure('Name', 'Absorption coefficient')
semilogx(general.freq, porous.alpha, '--k', 'LineWidth', 3); hold on;
semilogx(general.freq, alpha_PWA, 'b', 'LineWidth', 2); hold on;
semilogx(general.freq, alpha_q, 'r', 'LineWidth', 2); hold on;
xlabel('Frequency [Hz]');
ylabel('\alpha [-]');
grid on;
legend('Reference', 'PWA', 'q-term')
ylim([-0.4 1.4])



