% demo for image registration

close all 
clear variables
clc

[X,Y] = meshgrid(-2:.2:2, -2:.2:2);                                
Z = X .* exp(-X.^2 - Y.^2); 
inG.Y = Z;
[nx,ny] = size(inG.Y);

phi = [2;-2;-1];
sigma = 1e6;
[gx] = g_rigid2D([],phi,[],inG);
y = gx +  (1/sqrt(sigma))*randn(size(gx));

hf = figure;
subplot(2,2,1),imagesc(inG.Y),colormap(bone),colorbar,axis equal, axis tight,title('original image')
subplot(2,2,2),imagesc(reshape(y,nx,ny)),colormap(bone),colorbar,axis equal, axis tight,title('disaligned image')

g_fname = @g_rigid2D;
dim.n = 0;
dim.n_theta = 0;
dim.n_phi = 3;
options.inG = inG;
[posterior,out] = VBA_NLStateSpaceModel(...
    y,...
    [],...
    [],...
    g_fname,...
    dim,...
    options);
hg = out.suffStat.gx;

figure(hf)
subplot(2,2,3),imagesc(reshape(hg,nx,ny)),colormap(bone),colorbar,axis equal, axis tight,title('realigned image')
subplot(2,2,4),imagesc(reshape(hg-y,nx,ny)),colormap(bone),colorbar,axis equal, axis tight,title('realigned - disaligned')