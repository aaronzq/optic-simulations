%%% Gerchberg-saxton phase retrival algorithm
clear 
% close all
addpath('utils\');

%%% parameters
N = 1024;
lambda = 500e-9; 
% dx1 = lambda/2.05;% asm
dx1 = 5e-6; %fresnel
L1 = N*dx1;
focal = 36e-3;
iters = 100;
radius_doe = 2.23e-3;
distance = linspace(35.5e-3,36.5e-3,5);
% distance = linspace(35.8e-3,37.2e-3,5);
% distance = 36e-3;
fresnel = L1^2/lambda/focal

% pupil function, the source
x1 = dx1*(ceil(-N/2):ceil(N/2)-1);
y1 = x1;
[X1,Y1] = meshgrid(x1,y1);

% define doe circular region
doe_region = (X1.^2+Y1.^2)<radius_doe^2;
rings=getRings(X1,Y1,radius_doe,length(distance));
% rings = rings(:,:,[2,3,1,4,5]);

source_amplitude = ones(N).*doe_region;
source_phase = zeros(N).*doe_region;
% figure; surf(source_amplitude); 

% image planes, targets
x0 = 0;     		% center
y0 = 0;     		% center
sigma = 1e-6; 			% beam waist
A = 1;      		% peak of the beam 
res = ((X1-x0).^2 + (Y1-y0).^2)./(2*sigma^2);
psf_intensity = A  * exp(-res);

target_amplitude = sqrt(psf_intensity);
target_amplitude = repmat(target_amplitude,[1,1,length(distance)]);
target_phase = zeros(N,N,length(distance));
% target_phase = (2*rand(N,N) - 1)*pi;
%%
error = [];
target_amplitude_update = target_amplitude;
target_phase_update = target_phase;
source_phase_update = source_phase;
figure;
ht = title(['Updated amplitude at iteration: ' num2str(1)]);
for i = 1:iters
    i
    for d = 1:length(distance)
        target_update = target_amplitude(:,:,d).*exp(1i*target_phase_update(:,:,d));
        
        source_update = fresnelTF2d(target_update, dx1, -distance(d), lambda);

        temp = angle(source_update);
        temp = rotation_rep(temp(round((N-1)/2+1),:));
        source_phase_update(rings(:,:,d)>0) = temp(rings(:,:,d)>0);
        
        source_update = source_amplitude.*exp(1i*source_phase_update);
        
        target_update = fresnelTF2d(source_update, dx1, distance(d), lambda); 

        target_phase_update(:,:,d) = angle(target_update);
        
%         if d == round(length(distance)/2)
%             imagesc(target_amplitude_update);
%             title(['Updated amplitude at iteration: ' num2str(i)]);
%             pause(0.001);    
%             error = [error; sqrt(sum((target_amplitude-target_amplitude_update).^2,'all'))]; 
%         end

    end
    PSFs = zeros(N,N,length(distance));
    for dd=1:length(distance)
        [temp, ~, ~] = fresnelTF2d(source_update, dx1, distance(dd), lambda);
        PSFs(:,:,dd) = abs(temp).^2;
        subplot(1,length(distance),dd); imagesc(PSFs(512-50:513+50,512-50:513+50,dd))
        pause(0.001); 
    end
end
%%
% figure; plot(error);
% % figure; subplot(2,2,1); imagesc(abs(source_update));
% % subplot(2,2,2); imagesc((angle(source_update))); 
% % subplot(2,2,3); imagesc(abs(target_update));
% % subplot(2,2,4); imagesc((angle(target_update)));
% 
figure; imagesc((angle(source_update))); 
figure; plot(unwrap(angle(source_update(512,:))));
% 
% lambda*focal^2/radius_doe^2*1e6
% 
% 
%%
distance = linspace(35.5e-3,36.5e-3,5)-0.2e-3;
PSFs = zeros(N,N,length(distance));
figure;
for d=1:length(distance)
[temp, ~, ~] = fresnelTF2d(source_update, dx1, distance(d), lambda);
PSFs(:,:,d) = abs(temp).^2;
subplot(1,length(distance),d); imagesc(PSFs(512-50:513+50,512-50:513+50,d))
end
PSFs = PSFs ./ max(PSFs(:));
figure; plot(distance,squeeze(PSFs(512,512,:)))
