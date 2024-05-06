function [f2, dx2, x2] = cameraSampling(f1, dx1, pixel_size, resolution)
% assuming dx1=dy1, Nx=Ny

[~, N] = size(f1);

x1space = dx1*[ceil(-N/2):ceil(N/2)-1];
[X1, Y1] = meshgrid(x1space, x1space);

dx2 = pixel_size;
x2space = dx2*[ceil(-resolution/2):ceil(resolution/2)-1];
[X2, Y2] = meshgrid(x2space, x2space);

% determine if camera framesize is smaller than simulated image
assert(x2space(end)<=x1space(end), sprintf("camera framesize: %f should be smaller than simulated image: %f. \n", x2space(end)*2, x1space(end)*2 ));

f2 = interp2(X1, Y1, f1, X2, Y2);
x2 = x2space;

end
