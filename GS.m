%%% Gerchberg-saxton phase retrival algorithm
clear all
% close all
addpath('utils\');

%%% parameters
N = 200;
lambda = 500e-9; 
% dx1 = lambda/2.05;% asm
dx1 = 2e-5; %fresnel
L1 = N*dx1;
distance = 100e-5;
iters = 100;

L1^2/lambda/distance

x1 = dx1*(ceil(-N/2):ceil(N/2)-1);
y1 = x1;
[X1,Y1] = meshgrid(x1,y1);
x0 = 0;     		% center
y0 = 0;     		% center
sigma = 1000e-3; 			% beam waist
A = 1;      		% peak of the beam 
res = ((X1-x0).^2 + (Y1-y0).^2)./(2*sigma^2);
input_intensity = A  * exp(-res);
figure; surf(input_intensity); shading interp

forward_func = @(f1) ifftshift(fft2(fftshift(f1)));
backward_func = @(f1) ifftshift(ifft2(fftshift(f1)));
% 
% forward_func = @(f1) fresnelTF2d(f1, dx1, distance, lambda);
% backward_func = @(f1) fresnelTF2d(f1, dx1, -distance, lambda);

% forward_func = @(f1) asm2d(f1, dx1, distance, lambda);
% backward_func = @(f1) asm2d(f1, dx1, -distance, lambda);

source_amplitude = input_intensity;

img2d = im2double(imresize(imread("cell.tif"),[N,N]));
img2d = img2d(:,:,1);
target_amplitude = sqrt(img2d);
target_phase = zeros(N);
% target_phase = (2*rand(N,N) - 1)*pi;
target = target_amplitude.*exp(1i*target_phase);

error = [];
target_phase_update = target_phase;
figure;
for i = 1:iters
    
    target_update = target_amplitude.*exp(1i*target_phase_update);

    source_update = backward_func(target_update);
    source_amplitude_update = abs(source_update);
%     max(source_amplitude_update(:))
    source_phase_update = angle(source_update);

    source_update = source_amplitude.*exp(1i*source_phase_update);

    target_update = forward_func(source_update);
    target_amplitude_update = abs(target_update);
%     max(target_amplitude_update(:))
    target_phase_update = angle(target_update);
    
    imagesc(target_amplitude_update);
    title(['Updated amplitude at iteration: ' num2str(i)]);
    colormap('gray')
    colorbar
    axis equal
    pause(0.0001);

    error = [error; sqrt(sum((target_amplitude-target_amplitude_update).^2,'all'))]; 
end

figure; plot(error);
figure; subplot(2,2,1); imagesc(abs(source_update));
subplot(2,2,2); imagesc(angle(source_update));
subplot(2,2,3); imagesc(abs(target_update));
subplot(2,2,4); imagesc(angle(target_update));

