function [f2, dx2, x2] = asm2d(f1, dx1, z, lambda)
% Wave propagation by the angular spectrum method (ASM)
% inputs: f1: complex amplitude wavefront
%         dx1: the pixel size of input wavefront
%         lambda: wavelength (meter)
% outputs: f2: wavefront after propagation
%          dx2: the pixel size of output wavefront
%          x2: the coordinate
% assuming dx1=dy1, Nx=Ny, square field
% written by Zhaoqiang Wang, tested on 2022a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 2*pi/lambda;
[~, N] = size(f1);

% the spatial frequency domain
du = 1/(N*dx1);
u = du*[ ceil(-N/2):ceil(N/2)-1 ];
[U,V] = meshgrid(u,u);

% the propagator for each angular components
mask = (U.^2 + V.^2)<(1/lambda)^2;

max(u)>1/lambda

H = mask.*exp(1i*real(2*pi/lambda*sqrt(1-(lambda*U).^2-(lambda*V).^2)*z));

% propagation
f2 = ifftshift(ifft2(fftshift(ifftshift(fft2(fftshift(f1))).*H)));
% f2 = ifft2(fftshift(fft2(f1)).*H);

dx2 = dx1;
x2 = dx2*[ ceil(-N/2):ceil(N/2)-1 ];
end