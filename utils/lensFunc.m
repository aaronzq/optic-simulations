function [f2, dx2, x2] = lensFunc(f1, dx1, focal, lambda, aperture_diameter)

[Ny, Nx] = size(f1);
x1 = dx1*[ ceil(-Nx/2):ceil(Nx/2)-1 ];
[X, Y] = meshgrid(x1, x1);
aperture=double((X.^2+Y.^2)<(aperture_diameter/2)^2); % Objective lens clear aperture 

fprintf("Sampling term dx/lambda %f should be smaller than the lens with F/# f/d %f. Adequately sampling: %d\n", dx1/lambda, focal/aperture_diameter, dx1/lambda<focal/aperture_diameter);


PF = exp(-1i*pi/lambda/focal*(X.^2+Y.^2)).*aperture;
f2 = f1 .* PF;
dx2 = dx1;
x2 = x1;

end