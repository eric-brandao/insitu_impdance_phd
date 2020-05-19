function [p_ext, u_extr, u_extz]=extended_field_mat(hs,z,r,f,dimm,rho0,c0,m,n,message)
%%% This function calculates the sound field for an extended reactive
%%% backed by a rigid wall. Sample is infinite. 
%%% Boundary conditions are working well and particle velocity beneath the
%%% sample is correct up to the 4th decimal. Precision increases if B is
%%% higher
B=10;
w=2*pi*f;
k0=w/c0;
R1=sqrt(r^2+(hs-z)^2); R2=sqrt(r^2+(hs+z)^2);

if length(f)>20
    h1 = waitbar(0,message);
end

for ms=1:length(f)
    
    if length(f)>20
        waitbar(ms/length(f))
    end
    
    %%%%%% Field above the sample
    if z>=0 %% Above the sample's surface
        %%% Pressure
        F=@(s)(((2.*exp(-k0(ms)*(sqrt(s.^2-1))*(hs+z))).*k0(ms).*s.*besselj(0,k0(ms).*s*r))./...
        (((sqrt(s.^2-1))+m(ms).*(sqrt(s.^2-n(ms).^2)).*tanh(k0(ms)*dimm*sqrt(s.^2-n(ms).^2)))));
        Ip=quadl(F,0,B);
%         Ip=integral(F,0,B,'RelTol',1e-8,'AbsTol',1e-13,'Waypoints',0:0.01:1)
        p_ext(ms)=((exp(-1i*k0(ms)*R1))/R1)-((exp(-1i*k0(ms)*R2))/R2)+Ip;

        %%% Particle velocity r
        G=@(s)(((2.*k0(ms).*(s.^2).*exp(-k0(ms)*(sqrt(s.^2-1))*(hs+z))).*besselj(1,k0(ms).*s*r))./...
        (((sqrt(s.^2-1))+m(ms).*(sqrt(s.^2-n(ms).^2)).*tanh(k0(ms)*dimm*sqrt(s.^2-n(ms).^2)))));
        Iu=quadl(G,0,B);
%         Iu=integral(G,0,B,'RelTol',1e-8,'AbsTol',1e-13,'Waypoints',0:0.0001:1);
        u_extr(ms)=(1./(rho0*c0)).*(((exp(-1i*k0(ms)*R1))/R1).*((1./(1i*k0(ms)*R1))+1).*((-r)/R1)-...
            ((exp(-1i*k0(ms)*R2))/R2).*((1./(1i*k0(ms)*R2))+1).*((-r)/R2)-(1/1i)*Iu);
%         clear Iu
        %%% Particle velocity z
        G=@(s)(((2.*(sqrt(s.^2-1)).*exp(-k0(ms)*(sqrt(s.^2-1))*(hs+z))).*k0(ms).*s.*besselj(0,k0(ms).*s*r))./...
        (((sqrt(s.^2-1))+m(ms).*(sqrt(s.^2-n(ms).^2)).*tanh(k0(ms)*dimm*sqrt(s.^2-n(ms).^2)))));
        Iu=quadl(G,0,B);
%         Iu=integral(G,0,B,'RelTol',1e-8,'AbsTol',1e-13,'Waypoints',0:0.0001:1);
%         u_extz(ms)=(1./(rho0*c0)).*(((exp(-1i*k0(ms)*R1))/R1).*((1./(1i*k0(ms)*R1))+1).*((hs-z)/R1)+...
%             ((exp(-1i*k0(ms)*R2))/R2).*((1./(1i*k0(ms)*R2))+1).*((hs+z)/R2)-(1/1i)*Iu);
        u_extz(ms)=(((exp(-1i*k0(ms)*R1))/R1).*((1./(1i*k0(ms)*R1))+1).*((hs-z)/R1)+...
            ((exp(-1i*k0(ms)*R2))/R2).*((1./(1i*k0(ms)*R2))+1).*((hs+z)/R2)-(1/1i)*Iu);
%         clear Iu
        
    else %% Below the sample's surface
        k1=n(ms)*k0*(ms);
        c1=2*pi*f(ms)/k1(ms);
        rho1=rho0./m(ms);
        %%% Pressure
        F=@(s)(((cosh(k0(ms)*(sqrt(s.^2-n(ms).^2)).*(z+dimm)))./(cosh(k0(ms)*(sqrt(s.^2-n(ms).^2)).*dimm))).*...
            (((2.*exp(-k0(ms)*(sqrt(s.^2-1))*hs)).*k0(ms).*s.*besselj(0,k0(ms).*s*r))./...
        (((sqrt(s.^2-1))+m(ms).*(sqrt(s.^2-n(ms).^2)).*tanh(k0(ms)*dimm*sqrt(s.^2-n(ms).^2))))));
        Ip=quadl(F,0,B);
%         Ip=integral(F,0,B,'RelTol',1e-8,'AbsTol',1e-13,'Waypoints',0:0.0001:1);
        p_ext(ms)=Ip;
%         clear Ip
        %%%%%%%% Particle velocity r
        G=@(s)(((cosh(k0(ms)*(sqrt(s.^2-n(ms).^2)).*(z+dimm)))./(cosh(k0(ms)*(sqrt(s.^2-n(ms).^2)).*dimm))).*...
            (((2.*k0(ms)*s.*exp(-k0(ms)*(sqrt(s.^2-1))*hs)).*k0(ms).*s.*besselj(1,k0(ms).*s*r))./...
        (((sqrt(s.^2-1))+m(ms).*(sqrt(s.^2-n(ms).^2)).*tanh(k0(ms)*dimm*sqrt(s.^2-n(ms).^2))))));

        Iu=quadl(G,0,B);                
        u_extr(ms)=-(1./(1i*w(ms)*rho1))*Iu;
%         Iu=integral(G,0,B,'RelTol',1e-8,'AbsTol',1e-13,'Waypoints',0:0.0001:1);

        %%%%%%%%% Particle velocity r by pressure Gradient
%         Delta=1e-5; 
%         G=@(s)(((cosh(k0(ms)*(sqrt(s.^2-n(ms).^2)).*(z+dimm)))./(cosh(k0(ms)*(sqrt(s.^2-n(ms).^2)).*dimm))).*...
%             (((2.*exp(-k0(ms)*(sqrt(s.^2-1))*hs)).*k0(ms).*s.*besselj(0,k0(ms).*s*(r+Delta)))./...
%         (((sqrt(s.^2-1))+m(ms).*(sqrt(s.^2-n(ms).^2)).*tanh(k0(ms)*dimm*sqrt(s.^2-n(ms).^2))))));
%         Ip=quadl(G,0,B);
%         u_extr(ms)=(1./(rho1*c1)).*(1./(1i*n(ms)*k0(ms))).*((Ip-p_ext(ms))/Delta);
%         u_extr(ms)=(1./(1i*w(ms)*rho1*Delta)).*(Ip-p_ext(ms));
%         clear Iu
        %%%%%%%% Particle velocity z
        G=@(s)(((k0(ms)*sinh(k0(ms)*(sqrt(s.^2-n(ms).^2)).*(z+dimm)))./(cosh(k0(ms)*(sqrt(s.^2-n(ms).^2)).*dimm))).*...
            (((2.*(sqrt(s.^2-n(ms).^2)).*exp(-k0(ms)*(sqrt(s.^2-1))*hs)).*k0(ms).*s.*besselj(0,k0(ms).*s*r))./...
        (((sqrt(s.^2-1))+m(ms).*(sqrt(s.^2-n(ms).^2)).*tanh(k0(ms)*dimm*sqrt(s.^2-n(ms).^2))))));
        Iu=quadl(G,0,B);
%         Iu=integral(G,0,B,'RelTol',1e-8,'AbsTol',1e-13,'Waypoints',0:0.0001:1);
        u_extz(ms)=(1./(1i*w(ms)*rho1)).*Iu;  

%         clear Iu

%       %%%%%%%%% Particle velocity r byPressure Gradient
%         Delta=1e-5; 
%         G=@(s)(((cosh(k0(ms)*(sqrt(s.^2-n(ms).^2)).*(z+Delta+dimm)))./(cosh(k0(ms)*(sqrt(s.^2-n(ms).^2)).*dimm))).*...
%             (((2.*exp(-k0(ms)*(sqrt(s.^2-1))*hs)).*k0(ms).*s.*besselj(0,k0(ms).*s*(r)))./...
%         (((sqrt(s.^2-1))+m(ms).*(sqrt(s.^2-n(ms).^2)).*tanh(k0(ms)*dimm*sqrt(s.^2-n(ms).^2))))));
%         Ip=quadl(G,0,B);
% %         u_extr(ms)=(1./(rho1*c1)).*(1./(1i*n(ms)*k0(ms))).*((Ip-p_ext(ms))/Delta);
%         u_extz(ms)=(1./(1i*w(ms)*rho1*Delta)).*(Ip-p_ext(ms));
    end    
end
if length(f)>20
    close(h1)
end




% %%%% Incident pressure and particle velocity above the sample
%     if z>=0 %% Above the sample's surface
%         p_i=(exp(-1i*k0*R1))/R1;
%         uir=(1./(rho0*c0)).*((exp(-1i*k0*R1))/R1).*((1./(1i*k0*R1))+1).*((-r)/R1);
%         uiz=(1./(rho0*c0)).*((exp(-1i*k0*R1))/R1).*((1./(1i*k0*R1))+1).*((hs-z)/R1);
%     end
% 
% if length(f)>20
%     h1 = waitbar(0,strcat('Calculating spectra for receiver: r=',num2str(r),', z=',num2str(z)));
% end
% % matlabpool(2);
% parfor ms=1:length(f)
% %     
%     if length(f)>20
%         waitbar(ms/length(f))
%     end
%     
%     %%%%%% Reflected Field above the sample
%     if z>=0 %% Above the sample's surface
% %         %%% Incident pressure above the sample
% %         pi(ms)=(exp(-1i*k0(ms)*R1))/R1;
% %         
%         %%% Reflected Pressure
%         F=@(s)(((2.*exp(-k0(ms)*(sqrt(s.^2-1))*(hs+z))).*k0(ms).*s.*besselj(0,k0(ms).*s*r))./...
%         (((sqrt(s.^2-1))+m(ms).*(sqrt(s.^2-n(ms).^2)).*tanh(k0(ms)*dimm*sqrt(s.^2-n(ms).^2)))));
% %         Ip=integral(F,0,B);
%         Ip=quadl(F,0,B);
% %         Ip=integral(F,0,B,'RelTol',1e-8,'AbsTol',1e-13);
%         pr(ms)=-((exp(-1i*k0(ms)*R2))/R2)+Ip;
% 
%         %%% Reflected particle velocity r
%         G=@(s)(((2.*k0(ms).*(s.^2).*exp(-k0(ms)*(sqrt(s.^2-1))*(hs+z))).*besselj(1,k0(ms).*s*r))./...
%         (((sqrt(s.^2-1))+m(ms).*(sqrt(s.^2-n(ms).^2)).*tanh(k0(ms)*dimm*sqrt(s.^2-n(ms).^2)))));
%         Iu=integral(G,0,B);
% %         Iu=quadl(G,0,B);
% %         Iu=integral(G,0,B,'RelTol',1e-8,'AbsTol',1e-13,'Waypoints',0:0.0001:1);
%         urr(ms)=(1./(rho0*c0)).*(-((exp(-1i*k0(ms)*R2))/R2).*((1./(1i*k0(ms)*R2))+1).*((-r)/R2)-...
%             (1/1i)*Iu);
% 
%         %%% Reflected particle velocity z
%         G=@(s)(((2.*(sqrt(s.^2-1)).*exp(-k0(ms)*(sqrt(s.^2-1))*(hs+z))).*k0(ms).*s.*...
%             besselj(0,k0(ms).*s*r))./(((sqrt(s.^2-1))+m(ms).*(sqrt(s.^2-n(ms).^2)).*...
%             tanh(k0(ms)*dimm*sqrt(s.^2-n(ms).^2)))));
%         Iu=integral(G,0,B);
% %         Iu=quadl(G,0,B);
% %         Iu=integral(G,0,B,'RelTol',1e-8,'AbsTol',1e-13,'Waypoints',0:0.0001:1);
%         urz(ms)=(1./(rho0*c0)).*(((exp(-1i*k0(ms)*R2))/R2).*((1./(1i*k0(ms)*R2))+1).*((hs+z)/R2)-...
%             (1/1i)*Iu);
% %         clear Iu
%         
%     else %%%%%% Transmited field (Incidente and Reflected inside the sample)
%         k1=n(ms)*k0*(ms);
%         c1=2*pi*f(ms)/k1(ms);
%         rho1=rho0./m(ms);
%         
%         %%% Incident pressure in the sample
%         F=@(s)(((k0(ms)*s.*besselj(0,k0(ms).*s*r))./(cosh(k0(ms)*(sqrt(s.^2-n(ms).^2)).*dimm))).*...
%             ((exp(-k0(ms)*((sqrt(s.^2-1))*hs-((sqrt(s.^2-n(ms)))*(z+dimm)))))./...
%         (((sqrt(s.^2-1))+m(ms).*(sqrt(s.^2-n(ms).^2)).*tanh(k0(ms)*dimm*sqrt(s.^2-n(ms).^2))))));
% %         F=@(s)(((k0(ms)*s.*besselj(0,k0(ms).*s*r).*tanh(k0(ms)*dimm*sqrt(s.^2-n(ms).^2)))./...
% %             (sinh(k0(ms)*(sqrt(s.^2-n(ms).^2)).*dimm))).*...
% %             ((exp(-k0(ms)*((sqrt(s.^2-1))*hs-((sqrt(s.^2-n(ms)))*(z+dimm)))))./...
% %         (((sqrt(s.^2-1))+m(ms).*(sqrt(s.^2-n(ms).^2)).*tanh(k0(ms)*dimm*sqrt(s.^2-n(ms).^2))))));
% %         Ip=integral(F,0,B);
%         Ip=quadl(F,0,B);
% %         Ip=integral(F,0,B,'RelTol',1e-8,'AbsTol',1e-13,'Waypoints',0:0.0001:1);
%         p_i(ms)=Ip;
% 
%         %%% Reflected pressure in the sample
%         F=@(s)(((k0(ms)*s.*besselj(0,k0(ms).*s*r))./(cosh(k0(ms)*(sqrt(s.^2-n(ms).^2)).*dimm))).*...
%             ((exp(-k0(ms)*((sqrt(s.^2-1))*hs+((sqrt(s.^2-n(ms)))*(z+dimm)))))./...
%         (((sqrt(s.^2-1))+m(ms).*(sqrt(s.^2-n(ms).^2)).*tanh(k0(ms)*dimm*sqrt(s.^2-n(ms).^2))))));
% %           F=@(s)(((k0(ms)*s.*besselj(0,k0(ms).*s*r).*tanh(k0(ms)*dimm*sqrt(s.^2-n(ms).^2)))./...
% %             (sinh(k0(ms)*(sqrt(s.^2-n(ms).^2)).*dimm))).*...
% %             ((exp(-k0(ms)*((sqrt(s.^2-1))*hs+((sqrt(s.^2-n(ms)))*(z+dimm)))))./...
% %         (((sqrt(s.^2-1))+m(ms).*(sqrt(s.^2-n(ms).^2)).*tanh(k0(ms)*dimm*sqrt(s.^2-n(ms).^2))))));  
% %         Ip=integral(F,0,B);
%         Ip=quadl(F,0,B);
% %         Ip=integral(F,0,B,'RelTol',1e-8,'AbsTol',1e-13,'Waypoints',0:0.0001:1);
%         pr(ms)=Ip;
%         
%         %%% Incident Particle velocity r
%         G=@(s)((((k0(ms)^2)*(s.^2).*besselj(1,k0(ms)*s*r))./(cosh(k0(ms)*(sqrt(s.^2-n(ms).^2))*dimm))).*...
%             ((exp(-k0(ms)*((sqrt(s.^2-1))*hs-(sqrt(s.^2-n(ms).^2))*(z+dimm))))./...
%             ((sqrt(s.^2-1))+m(ms).*(sqrt(s.^2-n(ms).^2)).*tanh(k0(ms)*dimm*sqrt(s.^2-n(ms).^2)))));
%         Iu=integral(G,0,B);
%         uir(ms)=-(1./(1i*w(ms)*rho1))*Iu;
%         
%         %%% Reflected Particle velocity r
%         G=@(s)((((k0(ms)^2)*(s.^2).*besselj(1,k0(ms)*s*r))./(cosh(k0(ms)*(sqrt(s.^2-n(ms).^2))*dimm))).*...
%             ((exp(-k0(ms)*((sqrt(s.^2-1))*hs+(sqrt(s.^2-n(ms).^2))*(z+dimm))))./...
%             ((sqrt(s.^2-1))+m(ms).*(sqrt(s.^2-n(ms).^2)).*tanh(k0(ms)*dimm*sqrt(s.^2-n(ms).^2)))));
%         Iu=integral(G,0,B);
%         urr(ms)=-(1./(1i*w(ms)*rho1))*Iu;
%         
%         %%% Incident Particle velocity z
%         G=@(s)((((k0(ms)^2)*s.*(sqrt(s.^2-n(ms).^2)).*besselj(0,k0(ms)*s*r))./...
%             (cosh(k0(ms)*(sqrt(s.^2-n(ms).^2))*dimm))).*...
%             ((exp(-k0(ms)*((sqrt(s.^2-1))*hs-(sqrt(s.^2-n(ms).^2))*(z+dimm))))./...
%             ((sqrt(s.^2-1))+m(ms).*(sqrt(s.^2-n(ms).^2)).*tanh(k0(ms)*dimm*sqrt(s.^2-n(ms).^2)))));
%         Iu=integral(G,0,B);
%         uiz(ms)=(1./(1i*w(ms)*rho1))*Iu;
%         
%         %%% Reflected Particle velocity z
%         G=@(s)((((k0(ms)^2)*s.*(sqrt(s.^2-n(ms).^2)).*besselj(0,k0(ms)*s*r))./...
%             (cosh(k0(ms)*(sqrt(s.^2-n(ms).^2))*dimm))).*...
%             ((exp(-k0(ms)*((sqrt(s.^2-1))*hs+(sqrt(s.^2-n(ms).^2))*(z+dimm))))./...
%             ((sqrt(s.^2-1))+m(ms).*(sqrt(s.^2-n(ms).^2)).*tanh(k0(ms)*dimm*sqrt(s.^2-n(ms).^2)))));
%         Iu=integral(G,0,B);
%         urz(ms)=-(1./(1i*w(ms)*rho1))*Iu;
%     end    
% end
% if length(f)>20
%     close(h1)
% end




