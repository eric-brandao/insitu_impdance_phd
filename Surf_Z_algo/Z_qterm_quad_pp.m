function Z_qterm_quad=Z_qterm_quad_pp(f,c0,hs,z1,z2,r,ZIS,Hw)
% Calculates surface impedance by q-term interation process

%Limit the frequency range to speed up the calculation
Range=1:1:length(f);       % Every 20th frequency points is taken
Range(f(Range)<50)=[];     % Frequencies below 100 Hz are not calculated
Range(f(Range)>12000)=[];   % Frequencies above 10000 Hz are not calculated

%% Geometry
R1=sqrt((r^2)+((hs-z1)^2));
R1l=sqrt((r^2)+((hs+z1)^2));

R2=sqrt((r^2)+((hs-z2)^2));
R2l=sqrt((r^2)+((hs+z2)^2));

Zf=ZIS;         % Input impedance

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
                intp1=@(q)((exp((-q.*k)/Zg)).*...
    ((exp(-1i*k.*(sqrt((r^2)+(hs+z1-1i*q).^2))))./...
    (sqrt((r^2)+(hs+z1-1i*q).^2))));

                p1=(exp(-1i*k*R1)/R1)+(exp(-1i*k*R1l)/R1l)-...
                    (2*k/Zg)*(quad(intp1,0,20));
                
                intp2=@(q)((exp((-q.*k)/Zg)).*...
    ((exp(-1i*k.*(sqrt((r^2)+(hs+z2-1i*q).^2))))./...
    (sqrt((r^2)+(hs+z2-1i*q).^2))));

                p2=(exp(-1i*k*R2)/R2)+(exp(-1i*k*R2l)/R2l)-...
                    (2*k/Zg)*(quad(intp2,0,20));
                
                fZg=Hw(b)-p1./p2;
            end
            if abs(fZg)<0.000001;
                break
            end

    end
    %Zfinn(b,1)=Zg1;
    Z_qterm_quad(b,1)=Zg1;
end
close(hqterm)
clear D Zf Zg Zg1 Zg2 fZg fZg1 fZg2 a b ans

%R(Range)=(ZcalFINN(Range)-1)./(ZcalFINN(Range)+1);
