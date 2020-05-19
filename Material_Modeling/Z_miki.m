function [alpha, Zs, Zp, kp]=Z_miki(freq,rho0,c0,di,resist, theta)
    X = freq/resist;
    w = 2*pi*freq;
    k0 = w/c0;
    % Characteristic
    Zp = rho0 * c0 * ((1+0.07*X.^-0.632)-...
        1i*(0.107*X.^-0.632));
    kp = -1i*k0.*((0.16*X.^-0.618)+...
        1i*(1+0.109*X.^-0.618));
    % Refraction
    n_index = kp/k0;
    theta_t = asin(sin(theta./n_index));
    kzp = kp.*cos(theta_t);
    Zs = -1i*Zp.*(kp./kzp).*cot(kzp*di);
    Vp = (Zs*cos(theta) - rho0*c0)./(Zs*cos(theta) + rho0*c0);
    alpha = 1 - (abs(Vp)).^2;
end
