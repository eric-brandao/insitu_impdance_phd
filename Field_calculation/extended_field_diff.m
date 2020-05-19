function pres_rec = extended_field_diff(s_coord, r_coord,...
    freq, di, rho0, c0, Zp, kp, message)
%%% This function calculates the sound field for an extended reactive
%%% backed by a rigid wall under diffuse field incidence. Sample is infinite. 
%%%% material
% cp = 2*pi*f./kp;
w=2*pi*freq;
k0=w/c0;
rhop = (Zp.*kp)./w;
m = rho0./rhop;
n = kp./k0;

B=10;
hs = s_coord(:,1);

% if length(freq)>20
%     h1 = waitbar(0,message);
% end

%%%% Initialize pressure matrix
pres_rec = zeros(length(r_coord(:,1)), length(freq));

%%% Initialize receptor loop
for jrec = 1:length(r_coord(:,1))
    %%%% Message on command window
    disp(strcat('Calculating sound pressure for receiver ', num2str(jrec)));
    %%%% Needed distances r and zr
    r = sqrt((s_coord(:,1) - r_coord(jrec,1)).^2 +...
        (s_coord(:,2) - r_coord(jrec,2)).^2);
    zr = r_coord(jrec,3);
    r1=sqrt(r.^2+(hs-zr).^2); 
    r2=sqrt(r.^2+(hs+zr).^2);
    %%%% Loop over frequencies
    rng(0) %% set seed for randomization
%     q = zeros(length(freq), length(s_coord(:,1)));
    h1 = waitbar(0,message);
    for ms=1:length(freq)
        %%% waitbar
        if length(freq)>20
            waitbar(ms/length(freq))
        end
        %%% random amps
        amp = 0.0002 + (20-0.0002)*rand(length(s_coord(:,1)),1);
        phase = 2*pi*rand(length(s_coord(:,1)),1);
        q = amp.*exp(1i*phase);
        %%% Integrand
%         F=@(s)(sum(q.*(2.*exp(-k0(ms)*(sqrt(s.^2-1)).*(hs+zr))).*...
%             k0(ms).*s.*besselj(0,k0(ms).*s.*r))./...
%         (((sqrt(s.^2-1))+m(ms).*(sqrt(s.^2-n(ms).^2)).*...
%         tanh(k0(ms)*di*sqrt(s.^2-n(ms).^2)))));
        
        F=@(s)(sum(q.*(2.*exp(-k0(ms)*(sqrt(s.^2-1)).*(hs+zr))).*...
            k0(ms).*s.*besselj(0,k0(ms).*s.*r))./...
        (((sqrt(s.^2-1))+m(ms).*(sqrt(s.^2-n(ms).^2)).*...
        tanh(k0(ms)*di*sqrt(s.^2-n(ms).^2)))));
        %%%% integration
%         Ip=quadl(F,0,B);
        Ip=integral(F,0,B,'RelTol',1e-8,'AbsTol',1e-13,'Waypoints',0:0.0001:1);
        %%% sound pressure
        pres_rec(jrec, ms)=sum(((exp(-1i*k0(ms)*r1))./r1)-((exp(-1i*k0(ms)*r2))./r2))+Ip;
    end %%% end of freq loop
    close(h1)
%     if length(freq)>20
%         close(h1)
%     end
end




