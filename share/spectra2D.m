function [e_x, e_y, km] = spectra2D(dat, Lx)
% [e_x, e_y, km] = spectra2D(dat, Lx) calculate the 2D spectrum of dat.
%   Lx is the size of the domain (assuming Ly = Lx). e_x is the
%   x-coordinate wave number, e_y is the power, km is the mean wave number.

    nsize = size(dat);
    nnx = nsize(1);
    nny = nsize(2);
    fw = fftshift(fft2(dat));
    pw = abs(fw).^2;
    ew = zeros(nnx,1);
    km = 0;
    kxc=floor(nnx/2)+1;
    kyc=floor(nny/2)+1;
    for ky=1:nny
        for kx=1:nnx
            ki = sqrt((ky-kyc).^2+(kx-kxc).^2);
            k = round(ki)+1;
            ew(k) = ew(k) + pw(kx,ky);
            km = km + pw(kx,ky).*ki;
        end
    end
    kx = 1:1:nnx/2;
    sca = 2*pi/Lx;
    e_x = sca.*(kx-1);
    e_y = ew(1:nnx/2);
    km = km./sum(pw(:)).*sca;
end