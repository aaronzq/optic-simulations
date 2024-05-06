function [f2, dx2, x2] = fresnelIR2d(f1, dx1, z, lambda)
% assuming dx1=dy1, Nx=Ny

k = 2*pi/lambda;
[~, N] = size(f1);

% determine if it oversamples
fprintf("dx: %f; lambda*z/L: %f; Adequately sampling: %d\n", dx1, lambda*z/N/dx1, dx1<lambda*z/N/dx1);

x1 = dx1 * [ceil(-N/2):ceil(N/2)-1];
[X, Y] = meshgrid(x1, x1);

h = 1/(1j*lambda*z)*exp(1j*k/(2*z)*(X.^2+Y.^2));
H = fft2(fftshift(h))*dx1*dx1;
U1=fft2(fftshift(f1));
U2=H.*U1;
f2=ifftshift(ifft2(U2));

dx2 = dx1;
x2 = dx2*[ ceil(-N/2):ceil(N/2)-1 ];

end