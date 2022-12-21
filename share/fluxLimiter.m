function res = fluxLimiter(d1,d2,d3)
    % flux limiter
    % This function is used to calculate buoyancy flux from w and t,
    % matching the method used in the NCAR LES
    % ---- Cees's kappa=1/3 scheme

    r = (d1-d2+1.e-100)./(d2-d3-1.e-100);
    res = (d2-d3).*max(0., min(r, min(1./6+1./3.*r, 1.)));

