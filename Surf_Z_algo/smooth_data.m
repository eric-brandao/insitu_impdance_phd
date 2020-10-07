function Hw_sm=smooth_data(f,nfft,FS,Hw)
%%% This function smooths data obtained in the frequency dommain;
%%% Microflown algorithm

G=1-0.4/(1.5e4)*f;
for a=1:length(f) 
    b=round([a-G(a)*a*nfft/FS; a+G(a)*a*nfft/FS]);
    try
        Hw_sm(a,:)=mean(Hw((b(1):b(2)),1));
    catch
        Hw_sm(a,:)=Hw(a,1);
    end
end