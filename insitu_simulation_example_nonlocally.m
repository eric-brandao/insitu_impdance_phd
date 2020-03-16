clc; clearvars; close all
%%%% A minimalistic example on the simulation of in situ impedance
%%%% measurement.

%% Set some general configurations
general.freq = 100:20:4000; %%% frequency vector
general.c0 = 343; %%% sound speed
general.rho0 = 1.21; %%% air density
general.k0 = 2*pi*general.freq/general.c0;
%% Set sound source(s) - Define coordinates
sources(1).coord = [0, 0, 0.3]; 
sources(2).coord = [1, 0, 1.0]; 

%% Set receivers - Define coordinates
receivers(1).coord = [0, 0, 0.01];
receivers(2).coord = [0, 0, 0.02];
%% Porous material simulation
porous.resist = 9209; %%% Flow resistivity
porous.porosity = 0.99; %%% porosity
porous.a_inf = 1.00; %%% tortuosity
porous.Lam = 300e-6; %%% Vischous characteristic length
porous.Lam_l = 600e-6; %%% Thermal characteristic length
porous.thickness = 40/1000; %%% Sample thichness (over rigid wall)
porous.theta = deg2rad(0); %%% Angle of incidence

[porous.Zs, porous.Zp, porous.kp, porous.rhop]=pmaterial_jcal(2*pi*general.freq,general.rho0,...
    porous.thickness*1000,porous.resist,porous.porosity,...
    porous.a_inf,porous.Lam,porous.Lam_l);

porous.m = general.rho0./porous.rhop; % Density index
porous.n = porous.kp./general.k0; % refraction index

% Reflection and absorption coefficients
porous.Vp = (porous.Zs*cos(porous.theta)-general.c0*general.rho0)./...
    (porous.Zs*cos(porous.theta)+general.c0*general.rho0);
porous.alpha = 1 - (abs(porous.Vp)).^2;

%% sound field simulation (We first simulate the acoustic field - In another stage we retrive the surface impedance)
for js = 1:length(sources)
    hs = sources(js).coord(3); %% source height
    for jrec = 1:length(receivers) 
        hr = receivers(jrec).coord(3); %% receiver height
        r_vec = sources(js).coord - receivers(jrec).coord;
        r = sqrt(r_vec(1)^2 + r_vec(2)^2); %% horizontal distance between source and receiver
        message = strcat('calculating for source: ', num2str(js), ' and receiver: ', num2str(jrec));
        [sources(js).receivers(jrec).p, sources(js).receivers(jrec).ur, sources(js).receivers(jrec).uz] = ...
                extended_field_mat(hs,hr,r,general.freq,porous.thickness,...
                general.rho0,general.c0, porous.m, porous.n, message);
    end
end
clear jf k0 hq hr hs r r_vec js jrec
%% Now, let us recover the surface impedance (a q-term example).
%%% Obs: In the framework of this example, we need to choose a given
%%% source-receiver pair and then do the calculations. For instance for
%%% source(1) and receiver(1), we'd have
measurement.Zm = (sources(1).receivers(1).p./sources(1).receivers(1).uz)/(general.rho0 * general.c0);
hs = sources(1).coord(3); %% source height
hr = receivers(1).coord(3); %% receiver height
r_vec = sources(1).coord - receivers(1).coord;
r = sqrt(r_vec(1)^2 + r_vec(2)^2); %% horizontal distance between source and receiver
%%%% First, assuming a specular reflection of spherical wave
[~,alpha_PWA, Z_PWA]=Z_PWA_oblique(2*pi*general.freq/general.c0,...
    hs,hr,r,measurement.Zm);
%%% Now the q-term model
Zq = Z_qterm_quad_2(general.freq,general.c0,hs,...
    hr,r,Z_PWA,measurement.Zm);
alpha_q = 1 - (abs((Zq - 1)./(Zq + 1))).^2;
clear hr hs r r_vec
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

%% plots spk

% figure('Name', 'pu')
% subplot(2,2,1)
% semilogx(general.freq, 20*log10(abs(measurement.p)), '--k', 'LineWidth', 3); hold on;
% semilogx(general.freq, 20*log10(abs(sources(1).receivers(1).p)), 'b', 'LineWidth', 2); grid on;
% 
% subplot(2,2,2)
% semilogx(general.freq, angle(measurement.p), '--k', 'LineWidth', 3); hold on;
% semilogx(general.freq, angle(sources(1).receivers(1).p), 'b', 'LineWidth', 2); grid on;
% 
% subplot(2,2,3)
% semilogx(general.freq, 20*log10(abs(measurement.uz)), '--k', 'LineWidth', 3); hold on;
% semilogx(general.freq, 20*log10(abs(general.rho0*general.c0*sources(1).receivers(1).uz)), 'b', 'LineWidth', 2); grid on;
% 
% subplot(2,2,4)
% semilogx(general.freq, angle(measurement.uz), '--k', 'LineWidth', 3); hold on;
% semilogx(general.freq, angle(sources(1).receivers(1).uz), 'b', 'LineWidth', 2); grid on;
% xlabel('Frequency [Hz]');
% ylabel('\alpha [-]');
% 
% legend('loc', 'non loc')