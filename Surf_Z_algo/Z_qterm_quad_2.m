function Z_qterm_quad=Z_qterm_quad_2(f,c0,hs,z,r,ZIS,Zr)
% Calculates surface impedance by q-term interation process

%%%%%%%%%%%%%%%% Input list %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f --> frequency vector
% c0 --> Value for the sound speed
% hs --> Source height: vertical distance from the monopole source to the 
%        surface of the sample
% z --> Sensor height: same definition as hs
% r --> Source to sensor horizontal separation (if they are on the same
%       line r=0, which means normal incidence)
% ZIS --> Fist estimative of surface impedance; vector of the same length
%         of f; It comes from the Image source model (Z_PWA_oblique)
% Zr --> The impedance measured by the PU sensor (e.g 1 cm above the
%        surface of the sample); vector of the same length
%        of f

%%%%%%%%%%%%%%%% Output list %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Z_qterm_quad --> Surface impedance calculated with the q-term iterative
%                  method

%Limit the frequency range to speed up the calculation
Range=1:1:length(f);       % Every 20th frequency points is taken
Range(f(Range)<100)=[];     % Frequencies below 100 Hz are not calculated
Range(f(Range)>12000)=[];   % Frequencies above 10000 Hz are not calculated

%% Geometry
R1=sqrt((r^2)+((hs-z)^2));
R2=sqrt((r^2)+((hs+z)^2));

Zf=ZIS;         % First Guess - It comes from a simpler model of the sound field

% Variables:
% Here          Name in Finn's paper
% D         :   D'
% Zf(b)     :   Zm (measured impedance at the sensor location)
% Zg        :   Z' or Zk+1 (next surface impedance guess)
% Zg1       :   Zk (guessed surface impedance)
% Zg2       :   Zk-1 (previous guessed surface impedance)
% fZg       :   f(Z') (difference between measured Zm and calucation with Zg)
% fZg1      :   f(Z') (difference between measured Zm and calucation with Zg1)
% fZg2      :   f(Z') (difference between measured Zm and calucation with Zg2)

Z_qterm_quad(:,1)=zeros(length(f),1);
hqterm = waitbar(0,'Calculating Z by q-term Model...');
for b=Range;
    waitbar((b)/max(Range),hqterm)
    k =2*pi*f(b)/c0; %Wavenumber
    for a=1:40;
            if a==1         %first iteration the guessed impedance is the measured impedance
                Zg=Zf(b);   Zg1=1;      fZg=1;
            elseif a==2;    %second iteration the guessed impedance is first calculated imp
                            Zg1=Zg;     fZg1=fZg;   Zg=Zg-fZg;
            else            %all other iterations the secant method is used to guess the impedance
                if isfinite(Zg) && Zg~=Zg1 && fZg1~=fZg % To avoid errors
                            Zg2=Zg1;    fZg2=fZg1;
                            Zg1=Zg;     fZg1=fZg;

                    Zg=Zg1-(Zg1-Zg2)/(fZg1-fZg2)*fZg1;
                end
            end

            if isfinite(Zg) && Zg~=Zg1 && abs(fZg)>0.000001; % To avoid errors
                % Integrand of particle velocity  
                intu=@(q)((exp((-q.*k)/Zg)).*...
    ((exp(-1i*k.*(sqrt((r^2)+(hs+z-1i*q).^2))))./...
    (sqrt((r^2)+(hs+z-1i*q).^2))).*...
    ((hs+z-1i*q)./sqrt(r^2+(hs+z-1i*q).^2)).*...
    ((1./(1i*k*(sqrt((r^2)+(hs+z-1i*q).^2))))+1));                             
                
                % Actual particle velocity calculation (uses adaptative
                % quadrature)
                D=((exp(-1i*k*R1)/R1)*(1/(1i*k*R1)+1))*((hs-z)/R1)-...
                    ((exp(-1i*k*R2)/R2)*(1/(1i*k*R2)+1))*((hs+z)/R2)+...
                    (2*k/Zg)*quad(intu,0,20);
                
                % Integrand of sound pressure
                intp=@(q)((exp((-q.*k)/Zg)).*...
    ((exp(-1i*k.*(sqrt((r^2)+(hs+z-1i*q).^2))))./...
    (sqrt((r^2)+(hs+z-1i*q).^2))));
                
                % Difference between the actual measured impedance (Zr) and
                % the caculated impedance at sensor´s position at each
                % iteration
                fZg=Zr(b)-...
                    (((exp(-1i*k*R1)/R1)+(exp(-1i*k*R2)/R2)-...
                    (2*k/Zg)*(quad(intp,0,20)))./D);
                clc;
            end
            % Check convergence and stop if YES
            if abs(fZg)<0.000001;
                break
            end
    end
    Z_qterm_quad(b,1)=Zg1;
end
close(hqterm)
clear D Zf Zg Zg1 Zg2 fZg fZg1 fZg2 a b ans

