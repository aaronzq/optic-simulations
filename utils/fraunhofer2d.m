function [f2, dx2, x2] = fraunhofer2d(f1, dx1, z, lambda)
% assuming dx1=dy1, Nx=Ny

k = 2*pi/lambda;
[~, N] = size(f1);

% determine if it oversamples
fprintf("dx: %f should be larger than lambda*z/L: %f; Adequately sampling: %d\n", dx1, lambda*z/N/dx1, dx1>lambda*z/N/dx1);

dx2 = lambda*z/dx1/N;
x2 = dx2*[ ceil(-N/2):ceil(N/2)-1 ];
[X2, Y2] = meshgrid(x2, x2);

c = exp(1j*k*z) * exp(1j*k/2/z*(X2.^2+Y2.^2)) / (1j*lambda*z);
f2 = c .* ifftshift(fft2(fftshift(f1))) * dx1^2;

end