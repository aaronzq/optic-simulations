clear all
%% Simulate talbot effect
addpath('./utils/');

%%% parameters
N = 200;
Nz = 500;
lambda = 500e-9; 
dx1 = 1.25e-6;
L1 = N*dx1;
grating_period = 125*lambda;
grating_width  = 5*lambda;

%%% propagation distance, could be defined by Talbot length
% prop_range = 2*grating_period^2/lambda; %%% approximate when lambda << grating period
prop_range = lambda/(1-sqrt(1-lambda^2/grating_period^2));
z = linspace(0,prop_range,Nz);

%%% define the grating
grating = zeros(N);
for i = 1:max(round(grating_width/dx1),1)
    grating(round(grating_period/2/dx1)+i-1:round(grating_period/dx1):N,:) = 1;
end


%%% Propagtion
x1 = dx1*[ceil(-N/2):ceil(N/2)-1];
[X1,Y1] = meshgrid(x1,x1);
Ph1 = zeros(N);
E1 = ones(N).*exp(1i*Ph1); % plane wave

cross = zeros(Nz,N);
for ii =1:Nz
    [E2, dx2, x2] = asm2d(E1.*grating, dx1, z(ii), lambda);
%     [E2, dx2, x2] = fresnelTF2d(E1.*grating, dx1, z(ii), lambda);
    cross(ii,:)=E2(:,round(N/2));
end

f1 = figure('Position',[524.3333333333333,429,1660,462]); 
% figure;
ax1=axes('Position', [-0.3 0.1 1 0.8]); imagesc(ax1,grating); axis equal; colormap('gray'); hold on
line([round(N/2),round(N/2)], [0,N], 'Color','red','LineStyle','--'); hold off
ax1.DataAspectRatio = [100 100 1];
ax1.PlotBoxAspectRatio = [1 1 1];
ax2=axes('Position', [0.15 0.1 1 0.8]); imagesc(ax2,imrotate(abs(cross).^2,90)); colorbar; 
ax2.DataAspectRatio = [100 100 1];
ax2.PlotBoxAspectRatio = [1 1 1];
