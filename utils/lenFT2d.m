function [f2, dx2, x2] = lenFT2d(f1, dx1, flens, lambda)
% assuming dx1=dy1, Nx=Ny

k = 2*pi/lambda;
[~, N] = size(f1);

dx2 = lambda*flens/dx1/N;
x2 = dx2*[ ceil(-N/2):ceil(N/2)-1 ];
[X2, Y2] = meshgrid(x2, x2);


Koi = -1i*exp(1i*k*flens)/(lambda*flens);
f2 = Koi .* ifftshift(fft2(fftshift(f1))) * dx1^2;

end