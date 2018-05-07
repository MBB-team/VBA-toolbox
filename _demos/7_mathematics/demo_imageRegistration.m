% demo for image registration

close all 
clear variables

% 1: create original 2D image
[X,Y] = meshgrid(-2:.2:2, -2:.2:2);                                
Z = X .* exp(-X.^2 - Y.^2); 
inG.Y = Z;
[nx,ny] = size(inG.Y);

% 2: dis-align original image with rotation & translation
phi = [2;-2;-1];
sigma = 1e6;
[gx] = g_rigid2D([],phi,[],inG);
y = gx +  (1/sqrt(sigma))*randn(size(gx));

% 3: display original and disaligned images
hf = figure;
subplot(2,2,1),imagesc(inG.Y),colormap(bone),colorbar,axis equal, axis tight,title('original image')
subplot(2,2,2),imagesc(reshape(y,nx,ny)),colormap(bone),colorbar,axis equal, axis tight,title('disaligned image (target)')

% 4: undo dis-alignement:
%   - data consist of original image
%   - alignement parameters are not known
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

% 5: display realigned and residuals images
figure(hf)
subplot(2,2,3),imagesc(reshape(hg,nx,ny)),colormap(bone),colorbar,axis equal, axis tight,title('realigned image')
subplot(2,2,4),imagesc(reshape(hg-y,nx,ny)),colormap(bone),colorbar,axis equal, axis tight,title('realigned - disaligned')

displayResults(posterior,out,y,[],[],[],phi,[],sigma);