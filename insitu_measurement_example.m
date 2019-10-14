clc; clearvars; format long; close all
%%% This is a post processing program for pu
%% Constants
general.c0=343; 
general.rho0=1.21;
general.nfft=8192;

%% source / receiver
source.hs = 0.33; 
source.r=0;
receiver.hr=0.015;  

%% Load calibration file. Calculates FF transfer funtion and correction factor
load('Calibration-22-Sep-2008.mat'); clc;
cal.xt=data.recdata;
[cal.TF, general.f]=tfestimate(cal.xt(:,2),cal.xt(:,1),[],[],general.nfft,settings.FS); %TF
general.w=2*pi*general.f; general.k0=general.w/general.c0;
cal.TF_smoothed=smooth_data(general.f,general.nfft,settings.FS,cal.TF);
cal.CF=((1i*general.k0*(source.hs-receiver.hr))./((1i*general.k0*(source.hs-receiver.hr))+1))./...
    cal.TF_smoothed;

%% Load measurement file. Calculates FF transfer funtion and correction factor
load('Measurement-22-Sep-2008-flamex.mat'); clc;
meas.xt=data.recdata;
[meas.TF,general.f]=tfestimate(meas.xt(:,2),meas.xt(:,1),[],[],general.nfft,settings.FS); %TF
general.w=2*pi*general.f; general.k0=general.w/general.c0;
meas.Zm=meas.TF.*cal.CF;    
meas.Zm=smooth_data(general.f,general.nfft,settings.FS,meas.Zm);

clear settings curdir data eventdata folder handles hObject k0 path pathdir savename w 
%% process the impedance
[~,alpha_PWA,Z_PWA]=Z_PWA_oblique(general.k0, source.hs,...
    receiver.hr, source.r, meas.Zm);

Zq = Z_qterm_quad_2(general.f, general.c0,source.hs,...
    receiver.hr,source.r, Z_PWA , meas.Zm);
alpha_q = 1 - (abs((Zq - 1)./(Zq + 1))).^2;

%% plots
figure('Name', 'Absorption coefficient')
semilogx(general.f, alpha_PWA, 'b', 'LineWidth', 2); hold on;
semilogx(general.f, alpha_q, 'r', 'LineWidth', 2); hold on;
xlabel('Frequency [Hz]');
ylabel('\alpha [-]');
grid on;
legend('PWA', 'q-term')
ylim([-0.4 1.4]); xlim([100 10000])

