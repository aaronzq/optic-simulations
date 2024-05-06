function [f2, dx2, x2] = fresnelTF1d(f1, dx1, z, lambda)
% assuming dx1=dy1, Nx=Ny

k = 2*pi/lambda;
N = length(f1);

% determine if it oversamples for TF operand
fprintf("Spatial sampling dx: %f should be larger than" + ...
    " propagation coefficient lambda*z/L: %f; Satisfied? : %d\n", dx1, lambda*z/N/dx1, dx1>lambda*z/N/dx1);

du = 1/(N*dx1);
u = du*[ ceil(-N/2):ceil(N/2)-1 ];

H = exp(-1j*pi*lambda*z*(u.^2));
H = fftshift(H);
f2 = exp(1j*k*z)*ifftshift( ifft2( fft2(fftshift(f1)) .* H ) );

dx2 = dx1;
x2 = dx2*[ ceil(-N/2):ceil(N/2)-1 ];

end

