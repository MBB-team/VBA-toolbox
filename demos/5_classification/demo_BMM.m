
% clear all
close all


% path_sampling = [fileparts(mfilename('fullpath')),filesep,'..',filesep,'sampling'];
% addpath(path_sampling)

% class frequencies:
alpha = ones(3,1);
alpha = alpha./sum(alpha);

% class patterns:
lambda = [0.95 0.01 0.5
          0.5 0.95 0.01
          0.01 0.5 0.95];
      
% sample size:
n = 5e3;

% Generate n data samples from a binomial mixture model (BMM):
[y,labels] = generateBMM(n,alpha,lambda);

figure,imagesc(y)
colormap(bone)
set(gcf,'color',[1 1 1])
set(gca,'ytick',1:3)
title('binary data')
xlabel('data samples')
ylabel('data dimensions')

% classify data samples using VB inversion of BMM

% number of classes:
K = 3;
% set priors:
% options.c = 1e0*ones(2,1);
% options.b = 1e0*ones(K,1);
% options.TolFun = 1e-5;
options.algo = 'VB';
% VB inverion routine
[PI,out] = MixtureOfBinomials(y,K,options);

% display results
hf = figure('color',[1 1 1]);
h(1) = subplot(2,2,1,'parent',hf);
imagesc(lambda,'parent',h(1)),colorbar
title(h(1),'Simulated pattern')
h(2) = subplot(2,2,2,'parent',hf);
imagesc(out.Elambda,'parent',h(2)),colorbar
title(h(2),'Estimated pattern')
h(3) = subplot(2,2,3,'parent',hf);
imagesc(labels,'parent',h(3)),colorbar
title(h(3),'Simulated labels')

h(4) = subplot(2,2,4,'parent',hf);
imagesc(PI,'parent',h(4)),colorbar
title(h(4),'Estimated labels')



